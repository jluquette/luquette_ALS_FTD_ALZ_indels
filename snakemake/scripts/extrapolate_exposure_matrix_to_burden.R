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
        snakemake@input['exposures'],
        snakemake@input['burdens'],
        snakemake@output['csv']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("usage: extrapolate_exposure_matrix_to_burden.R expomat.csv mutburdens.csv out.csv")
}

expomat.file <- args[1]
mutburden.file <- args[2]
out.csv <- args[3]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))

mutburden <- fread(mutburden.file)

mutburden[, scaling.factor := ifelse(nsom == 0, 0, genome.burden/nsom)]
setkey(mutburden, sample)

# Exposure matrix format: rows are samples, columns are signatures
# Samples   SBS96A  SBS96B  SBS96C  SBS96D  SBS96E  SBS96F
# 11901106_Gy_C1  328 12  37  19  34  19
# 11901106_Gy_C2  271 27  25  37  44  46
# 1278-Oligo-5    6   23  0   3   5   3
# 1278-Oligo-7    7   26  6   15  14  4
expomat <- fread(expomat.file)
orig.colnames <- colnames(expomat)
print(head(expomat))
expomat <- as.matrix(expomat, rownames=1)
print(head(expomat))
expomat <- expomat * mutburden[rownames(expomat), scaling.factor]
print(head(expomat))

# reconstruct data table from matrix
new.dt <- data.table(rownames(expomat), expomat)
colnames(new.dt) <- orig.colnames

fwrite(new.dt, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
