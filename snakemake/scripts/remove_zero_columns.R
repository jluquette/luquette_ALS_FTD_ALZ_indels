#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['mat'],
        snakemake@output['mat']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}


args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
    stop("usage: remove_zero_columns.R input_matrix.txt output_matrix.txt")
}

infile <- args[1]
outfile <- args[2]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))

suppressMessages(library(data.table))

mat <- fread(infile)

zero.cols <- 1 + which(colSums(mat[,-1]) == 0)
print(zero.cols)

fwrite(mat[,-..zero.cols], file=outfile, sep='\t', quote=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
