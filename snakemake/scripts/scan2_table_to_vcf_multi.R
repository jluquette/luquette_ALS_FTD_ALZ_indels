#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['csv'],
        snakemake@params['qualtype'],
        snakemake@params['filter'],
        snakemake@output['vcf'],
        unlist(snakemake@params['samples'])
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
} else {
    cat("This script must be run within snakemake. It does not support command line usage.\n")
}

incsv <- snakemake@input[['csv']]
qualtype <- snakemake@params[['qualtype']]
filter <- tolower(snakemake@params[['filter']])
# vcfs and samples must have a 1-to-1 mapping
outvcfs <- snakemake@output[['vcfs']]
samples <- snakemake@params[['samples']]

cat("incsv is", incsv, '\n')
cat("qualtype is", qualtype, '\n')
cat('filter is', filter, '\n')
str(outvcfs)
str(samples)

if (length(outvcfs) != length(samples))
    stop("output$vcfs and params$samples must match in length and be 1-to-1.")

if (!(qualtype %in% c('A', 'AB', 'indel_A', 'indel_AB')))
    stop('qualtype must be one of A, AB, indel_A or indel_AB (case sensitive)')

if (!(filter %in% c('filtered', 'unfiltered')))
    stop(paste('argument 3 must be either "filtered" or "unfiltered", got', filter))

if (any(sapply(outvcfs, file.exists)))
    stop(paste('an output file already exists, please delete it first'))

suppressMessages(library(data.table))

muttab <- fread(incsv)

# outdated name scheme for mutation calls: A = VAF-based (pass), B = mutsig rescued (rescue)
# indel_ for indels; no prefix means SNV
mt <- ifelse(substr(qualtype, 1, 5) == 'indel', 'indel', 'snv')
add.rescue <- substr(qualtype, nchar(qualtype), nchar(qualtype)) == 'B'

muttab <- muttab[muttype == mt & (pass == TRUE | (add.rescue & rescue == TRUE))]

if (!is.null(samples))
    muttab <- muttab[sample %in% samples]

if (filter == 'filtered')
    muttab <- muttab[final.filter == FALSE]


for (i in 1:length(outvcfs)) {
    this.sample <- samples[i]
    outvcf <- outvcfs[i]
    cat(this.sample, '->', outvcf, '\n')

    # dummy header
    lines.to.write <- c("##fileformat=VCFv4.0", "##source=scan2_table_to_vcf.R", 
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    lines.to.write <- c(lines.to.write, paste(c("#CHROM", "POS", "ID", 
        "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", this.sample), 
        collapse = "\t"))

    this.muttab <- muttab[sample == this.sample]

    if (nrow(this.muttab) > 0)
        lines.to.write <- c(lines.to.write, paste(this.muttab$chr, this.muttab$pos, ".", this.muttab$refnt, this.muttab$altnt, ".", "PASS", ".", 'GT', '0/1', sep='\t'))

    writeLines(lines.to.write, con=outvcf)
}

if ('snakemake' %in% ls()) {
    sink()
}
