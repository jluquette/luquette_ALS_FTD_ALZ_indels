#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=16G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['gtf'],
        snakemake@params['tissue'],
        snakemake@input['gct'],
        snakemake@params['transcript_type'],
        snakemake@output['raw_bigwig'],
        snakemake@output['adj_bigwig'],
        snakemake@output['plot']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 7) {
    stop('usage: make_gene_length_bigwig.R gencode_gene_model.gtf tissue_without_spaces_or_parens tissue_median_tpm.gct transcript_type out.raw_width.bigwig out.adjusted_width.bigwig diagnostic_plot.pdf')
}

gene.model.file <- args[1]
tissue <- args[2]
tissue.median.tpm.file <- args[3]
transcript.type <- args[4]
out.raw.bigwig <- args[5]
out.adj.bigwig <- args[6]
out.pdf <- args[7]

tile.size=50

if (file.exists(out.raw.bigwig))
    stop(paste('output file', out.raw.bigwig, 'already exists, please delete it first'))
if (file.exists(out.adj.bigwig))
    stop(paste('output file', out.adj.bigwig, 'already exists, please delete it first'))
if (file.exists(out.pdf))
    stop(paste('output file', out.pdf, 'already exists, please delete it first'))


suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))


cat('Loading gene model:', gene.model.file, '\n')
# this version is the result of running GTEx's collapse script
gtex <- import.gff(gene.model.file)
# ENSEMBL IDs have underscores after them in GENCODE but not in GTEx
gtex$gene_id2 <- sapply(strsplit(gtex$gene_id,'_'),head,1)

cat('Annotating tissue median TPM expression levels:', tissue.median.tpm.file, '\n')
gtex.tpm <- fread(tissue.median.tpm.file, skip=2)
setkey(gtex.tpm, Name)
# making the entire pipeline safe for file names with spaces/parens is WAY too much work
colnames(gtex.tpm) <- gsub(' ', '_', colnames(gtex.tpm))
colnames(gtex.tpm) <- gsub('(', '_', colnames(gtex.tpm), fixed=TRUE)
colnames(gtex.tpm) <- gsub(')', '_', colnames(gtex.tpm), fixed=TRUE)
gtex.tpm <- gtex.tpm[gtex$gene_id2]
gtex$mean.expr <- rowMeans(gtex.tpm[,-(1:2)])

cat(paste0('Using tissue=', tissue, ' to compute max expression signal along genome\n'))
cat('Min. non-zero tissue-specific expression = ', min(gtex.tpm[[tissue]][gtex.tpm[[tissue]]>0], na.rm=T), '\n')
eps=1e-4
cat("Adding epsilon =", eps, "to expression levels to avoid log(0)\n")
gtex$tissue.expr <- gtex.tpm[[tissue]]+eps

# mean.expr is used to remove suspect transcripts, thus does not use tissue specific information
selected.genes <- !is.na(gtex$mean.expr) & apply(gtex.tpm[,-(1:2)],1,max) > 0

cat(sprintf('    %d / %d genes with non-NA expression and median expression > 0 in at least 1 of %d tissues\n',
    sum(selected.genes), length(selected.genes), ncol(gtex.tpm)-2))

gtex <- gtex[selected.genes,]
#e <- e[selected.genes]

# Only retain the specified tissue_type
# N.B. manual exploration revealed that about 400 genes of ~1 kb in length have very low expression,
# leading to a strong length-expression bias. These genes are primarily named ORxxx and KRTAPxxx and
# appear to represent unusual transcripts. These are protein_coding genes. We remove them from
# consideration here.
gtex <- gtex[gtex$transcript_type == transcript.type  & !(substr(gtex$gene_name,1,2)=='OR' | substr(gtex$gene_name,1,5)=='KRTAP'),]
# Only retain primary contigs
gtex <- gtex[seqnames(gtex) %in% seqlevels(gtex)[1:25],]
seqlevels(gtex) <- seqlevels(gtex)[1:25]
seqinfo(gtex) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
gtex$length <- width(gtex)  # just to get it into the mcols(gtex) data.frame for `lm`

gtex

# Use a linear model to adjust the width to account for gene length-expression correlation
m <- lm(log10(length) ~ log10(tissue.expr), data=mcols(gtex))
summary(m)
gtex$adj.length <- resid(m)
q <- function(score, nquantiles=10) findInterval(score, quantile(score, na.rm=T, probs=1:nquantiles/nquantiles), rightmost.closed=TRUE)+1
#gtex
p <- (ggplot(as.data.frame(mcols(gtex)), aes(x=tissue.expr, y=length)) + geom_point(size=1/10) + geom_smooth(method='lm', se=F) + scale_x_log10() + scale_y_log10() + labs(title='Before regression')) |
     (ggplot(as.data.frame(mcols(gtex)), aes(x=tissue.expr, y=adj.length)) + geom_point(size=1/10) + geom_smooth(method='lm', se=F) + scale_x_log10() + labs(title='After regression')) |
    (ggplot(as.data.frame(mcols(gtex)), aes(x=q(length), y=q(adj.length))) + geom_jitter(size=1/10) + labs(title='Effect on quantile binning')) & theme(aspect.ratio=1) #+ plot_annotation(title=paste0('Transcript type: ', transcript.type))
ggsave(p, filename=out.pdf)

cat('Tiling hg19 with', tile.size,  'bp bins\n')
bins <- tileGenome(seqlengths=seqlengths(gtex), tilewidth=tile.size, cut.last.tile.in.chrom=T)

cat(paste0('Mapping GTEx to non-overlapping ', tile.size, ' bp bins and taking max raw gene length across genes at each site\n'))
ols <- findOverlaps(bins, gtex)
#system.time(maxexpr <- sapply(split(e[to(ols)], from(ols)), max))
system.time(maxlen <- sapply(split(gtex$length[to(ols)], from(ols)), max))
map1 <- data.frame(from=as.integer(names(maxlen)), to=maxlen)

cat(paste0('Mapping GTEx to non-overlapping ', tile.size, ' bp bins and taking max adjusted gene length across genes at each site\n'))
system.time(maxadjlen <- sapply(split(gtex$adj.length[to(ols)], from(ols)), max))
map2 <- data.frame(from=as.integer(names(maxadjlen)), to=maxadjlen)

# Make a fresh tiling GRanges with no metadata and attach the expression values for each tile.
new.bins <- GRanges(seqnames=seqnames(bins), ranges=ranges(bins), seqinfo=seqinfo(bins))
new.bins$score <- NA
new.bins$score[map1$from] <- map1$to
new.bins <- new.bins[!is.na(new.bins$score),]

cat("Writing raw length bigwig ", out.raw.bigwig, '\n')
export.bw(new.bins, con=out.raw.bigwig)

new.bins <- GRanges(seqnames=seqnames(bins), ranges=ranges(bins), seqinfo=seqinfo(bins))
new.bins$score <- NA
new.bins$score[map2$from] <- map2$to
new.bins <- new.bins[!is.na(new.bins$score),]
cat("Writing adjusted length bigwig ", out.adj.bigwig, '\n')
export.bw(new.bins, con=out.adj.bigwig)

if ('snakemake' %in% ls()) {
    sink()
}
