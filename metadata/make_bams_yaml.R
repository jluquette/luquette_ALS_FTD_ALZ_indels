#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

immutable_meta.csv <- args[1]
out.yaml <- args[2]
if (length(args) != 2) 
    stop('usage: make_bams_yaml.R immutable_metadata.csv out.yaml')

if (file.exists(out.yaml)) {
    stop(paste0('output file ', out.yaml, ' already exists, please delete it first'))
}

x <- read.table(immutable_meta.csv, header=T, sep=',')

lines <- c()
for (d in unique(x$donor)) {
    lines <- c(lines, paste0("'", d, "':"))
    xd <- x[x$donor == d & x$amp != 'bulk',]
    lines <- c(lines, paste0("    single_cell:"))
    for (s in unique(xd$sample))
        lines <- c(lines, paste0('        ', s, ': bams/', s, '.bam'))
    xd <- x[x$donor == d & x$amp == 'bulk',]
    lines <- c(lines, paste0("    bulk:"))
    for (s in unique(xd$sample))
        lines <- c(lines, paste0('        ', s, ': bams/', s, '.bam'))
}

writeLines(text=lines, con=out.yaml)
