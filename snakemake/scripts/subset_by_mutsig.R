#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv_or_rda'],
        snakemake@input['csv_or_rda'],
        snakemake@params['filetype'],
        snakemake@params['colname'],
        snakemake@params[['channels']]
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat("Subsets either a permutation set (.rda file) or mutation set (.csv file) to specific signature channels.\n")
    cat("If an .rda file is supplied, it must contain a single object of class GenomicRanges or CompressedGRangesList.\n")
    cat('the object in in.rda must contain a column named `mutsig` that contains SBS96 channel info\n')
    stop('usage: subset_by_mutsig.R out.{rda or csv} in.{rda or csv} {rda|csv} colname channel1 [ channel2 ... channelN ]')
}


out.file <- args[1]
in.file <- args[2]
file.type <- tolower(args[3])
colname <- args[4]
channels <- args[-(1:4)]

if (file.exists(out.file))
    stop(paste('output file', out.file, 'already exists, please delete it first'))


if (file.type == 'csv') {
    suppressMessages(library(data.table))
    cat('CSV file - mutation table\n')
    mut.tab <- fread(in.file)
    mut.tab <- mut.tab[mut.tab[[colname]] %in% channels]
    fwrite(mut.tab, file=out.file)
} else if (file.type == 'rda') {
    suppressMessages(library(GenomicRanges))
    o.name <- load(in.file)
    original <- get(o.name)
    if ('GenomicRanges' %in% class(original) | 'CompressedGRangesList' %in% class(original)) {
        cat('detected permutation list\n')
        subsetted <- GenomicRanges::GRangesList(lapply(original, function(x) x[mcols(x)[[colname]] %in% channels,]))
        # Saves the object named o.name, not the object o.name itself
        assign(o.name, subsetted)
        save(list=o.name, file=out.file)
    } else {
        stop(paste('expected GenomicRanges or CompressGRangesList for object', o.name, 'in .rda file, but got', class(original)))
    }
} else {
    stop("file type must be either 'csv' or 'rda'")
}


if ('snakemake' %in% ls()) {
    sink()
}
