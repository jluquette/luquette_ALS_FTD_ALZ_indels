#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    if (!('prevburdens' %in% names(snakemake@input)))
        snakemake@input['prevburdens'] <- 'Notafile'

    commandArgs <- function(...) unlist(c(
        snakemake@input['mutburden'],
        snakemake@input['expomat'],
        snakemake@params['sig_to_subtract'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    cat("expomat.csv is the matrix of signature exposures. It *MUST* have been computed with Signature B present or else this correction will not succeed.\n")
    cat("Furthermore, signature B exposure should *NOT* be present in expomat.csv.\n")
    cat("UPDATE: this script probably does not work from the command line because signature_to_subtract must be missing to revert to legacy behavior.\n")
    stop("usage: correct_mutburden_mda.R mutburden.csv expomat.csv signature_to_subtract out.csv")
}

mutburden.file <- args[1]
expomat.file <- args[2]
sig.to.subtract <- args[3]
out.csv <- args[4]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

mutburden <- fread(mutburden.file)
# Rename the current genome.burden column
mutburden[, uncorrected.genome.burden := genome.burden]
mutburden$genome.burden <- NULL


expomat <- fread(expomat.file)

if (is.na(sig.to.subtract)) {
    # Old format: columns are samples, rows are signatures
    # Sig,1465-cortex_1-neuron_MDA_12,1465-cortex_1-neuron_MDA_18,1465-cortex_1-neuron_MDA_20
    # SBS1,0,1.03213469283357,7.92262828383771
    # SBS5,96.0582711821305,189.947280050293,265.6563961136
    # SBS16,128.01683687437,90.8979596350206,21.3788055929376
    # SBS19,0,73.943601661218,7.68544515576892
    # SBS30,121.625500453215,119.179520346402,65.9800619148647
    # SBS32,20.6888273171965,101.62274827524,157.669916713619
    # SBS40a,0,0,0
    expo.sums <- colSums(expomat[,-1])
    expo.sums <- data.table(sample=names(expo.sums), genome.burden=expo.sums)
} else {
    # New format: rows are samples, columns are signatures
    # Samples   SBS96A  SBS96B  SBS96C  SBS96D  SBS96E  SBS96F
    # 11901106_Gy_C1  328 12  37  19  34  19
    # 11901106_Gy_C2  271 27  25  37  44  46
    # 1278-Oligo-5    6   23  0   3   5   3
    # 1278-Oligo-7    7   26  6   15  14  4
    expo.sums <- as.matrix(expomat, rownames=1)
print(head(expo.sums))
print(-which(colnames(expo.sums) == sig.to.subtract))
    expo.sums <- rowSums(expo.sums[,-which(colnames(expo.sums) == sig.to.subtract)])
print(head(expo.sums))
    expo.sums <- data.table(sample=names(expo.sums), genome.burden=expo.sums)[sample %in% mutburden$sample]
print(head(expo.sums))
}


# Corrected genome burden now occupies the "genome.burden" column.
mutburden <- mutburden[expo.sums,,on=.(sample)]

fwrite(mutburden, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
