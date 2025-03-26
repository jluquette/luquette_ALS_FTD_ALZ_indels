library(scan2)
library(mclust)
library(Biostrings)

scan2.objects <- list.files(path = "scan2/summary_objects", pattern = ".rda$", full.names=TRUE)

# Loads all 575 objects
# Get all raw VAF-based calls. These include some outlier cells and
# MDA artifact calls.
all.passing <- rbindlist(lapply(scan2.objects, function(f) {
    print(f)
    load(f)
    passing(smry)  # don't use muttype='indel' to avoid 0 indel cases
}))[rescue == FALSE & muttype == 'indel']

meta <- fread("metadata/sample_metadata.csv")

# Join metadata like amplification type (MDA vs PTA), age, phenotype..
all.passing <- meta[all.passing,, on=.(sample)]


mb <- do.call(rbind, lapply(scan2.objects, function(f) {
    print(f)
    load(f)
    matrix(mutburden(smry), nrow=1, dimnames=list(name(smry), c('snv', 'indel')))
}))
mb <- cbind(meta[amp!='bulk'], mb[meta[amp != 'bulk']$sample,])


# The "pink" colored peaks corresponding to 2bp deletions in repeats
pink.peaks <- c(
    paste0('2:Del:R:', 0:1)
)

# The purple peaks corresponding to 2bp and 3bp deletions with local microhomology
purple.peaks <- c(
    '2:Del:M:1',
    paste0('3:Del:M:', 1:2)
)

mb <- all.passing[, .(n.indels=length(mutsig), n.pink=sum(mutsig %in% pink.peaks), n.purple=sum(mutsig %in% purple.peaks)),by=.(sample)][, c('ratio.pink', 'ratio.purple') := list(n.pink/n.indels, n.purple/n.indels)][mb,,on=.(sample)]
mb[is.na(n.indels), c('n.indels', 'n.pink', 'n.purple', 'ratio.pink', 'ratio.purple') := list(0, 0, 0, 0, 0)]



# From lots of manual inspection
# XXX: TODO: there are 2 types of outliers:
#   1. mutburden outliers, which may only indicate a failure of extrapolation
#      (which can happen with extreme allele balance or depth distributions),
#      but for which the actual mutation calls are fine.
#   2. call outliers, where the number of VAF-based calls are out of line
outliers <- c(
    # All PFC (i.e., non-DG) neurons from UMB5532 were indel HIGH outliers,
    # supporting a batch effect explanation.
    paste0('5532pfc-', c('Lp1C5', 'Lp2C5', 'Rp1C11', 'Rp1H7')),
    # Only 1 of the 5 DG neurons from UMB5532 was an indel HIGH outlier.
    '5532-1cp1E9',
    # One odd DG neuron from 5559. 2 other DG neurons were fairly high (like 2x),
    # others (PFC and DG) were comparable.
    '5559-dg1C11',
    # These 4 ALZ neurons shared thousands of PASS calls; a 5th neuron that appears to
    # be from the same batch (MA1647_p1E05) shares only 2.  Consulted with Mike
    # Miller, who said these neurons were suspected to be contaminated by Psomagen
    # and were thus excluded from his analysis.
    paste0('MA1647_p1', c('B','B','H','C'), '0', c(6,5,5,5)),
    # The next 4 cells have very unusual SNV signatures dominated by GCN>GTN. This is
    # so prevalent that it dwarfs MDA artifact signatures even by an order of magnitude
    # and introduces an entirely new signature into the de novo output.
    #
    # It isn't clear whether this signature is artifactual. Donor MADRC1456 is male,
    # so it could be useful to check chrX.
    #
    # These 4 cells also have high levels of ID83C, a very odd indel signature and
    # extremely high SNV passA counts.
    '4556_170602_P1F05',
    '4556_170602_P1F07',
    '1456_180403_Hp1D12',
    '1647-181116-Hp3C03'
)


# UPDATE: Feb. 25, 2025: the default mclustBIC function now returns a different
# clustering as optimal. Before, the VII,4 model had the best BIC, but this is
# now the second best BIC after VEE,6.  The difference in BIC is minimal, but
# there are some minor changes in classification.
#
# I assume this is due to differences in the randomness of cluster fitting.
#
# Top 3 models based on the BIC criterion: 
#    VEE,6    VII,4    VEI,7 
# 2129.430 2126.888 2122.291 
bic1 <- mclustBIC(mb[n.indels>=20 & !(sample %in% outliers),.(ratio.pink,ratio.purple)])
mod1.2 <- Mclust(mb[n.indels>=20 & !(sample %in% outliers),.(ratio.pink,ratio.purple)], x=bic1, modelNames='VII')
mb[, mclust.class2 := 0]
mb[n.indels>=20 & !(sample %in% outliers), mclust.class2 := mod1.2$classification]
mb[, mclust.uncert2 := 0]
mb[n.indels>=20 & !(sample %in% outliers), mclust.uncert2 := mod1.2$uncertainty]


ggplot(mb, aes(x=ratio.pink, y=ratio.purple, shape=factor(celltype), col=factor(mclust.class2))) + geom_point()
