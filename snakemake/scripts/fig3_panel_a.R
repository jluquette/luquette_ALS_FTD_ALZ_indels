#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input[1:2], snakemake@output[1:3]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
    stop("usage: plot_raw_indel.R neuron.rda oligo.rda out.svg out.pdf out.csv")
}

nfile <- args[1]
gfile <- args[2]
outsvg <- args[3]
outpdf <- args[4]
outcsv <- args[5]

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(data.table))
suppressMessages(library(svglite))


figwidth=7
figheight=4
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    layout(1:2)
    par(mar=c(1,4,3,1))
    n <- fread(nfile)
    g <- fread(gfile)
    plot.id83(x=1, spectrum=n$Spectrum, main='Neurons')
    plot.id83(x=1, spectrum=g$Spectrum, main='Oligo')
}

d <- data.table(MutType=names(n$MutType), Neurons=as.vector(n$Spectrum), Oligo=as.vector(g$Spectrum))
fwrite(d, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
