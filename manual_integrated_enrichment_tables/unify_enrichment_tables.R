#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    stop('this script is not updated for snakemake')
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@input['csv'],
        snakemake@params['tags'],
        snakemake@input['expomats']       # variable length list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat("WARNING: metadata is mined from the names of enrichment_table.csv. Expected format: <group_name>___<mutation_type>.<ignored...>.csv\n")
    stop("usage: unify_enrichment_tables.R out.csv group_metadata_table2.csv enrichment_table1.csv [ enrichment_table2.csv ... enrichment_tableN.csv ]")
}

out.csv <- args[1]
groupmeta.csv <- args[2]
enrichment.csvs <- args[-(1:2)]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))

gmeta <- fread(groupmeta.csv)


firstfield <- sapply(strsplit(basename(enrichment.csvs), split=".", fixed=TRUE), head, 1)
groupnames <- sapply(strsplit(firstfield, split="___"), `[`, 1)
mutnames <- sapply(strsplit(firstfield, split="___"), `[`, 2)


tab <- rbindlist(lapply(seq_along(enrichment.csvs), function(i) {
    fname <- enrichment.csvs[i]
    print(fname)
    x <- fread(fname)
    # some enrichment data types have a "group" metadata variable already. rename it.
    colnames(x)[colnames(x)=='group'] <- 'enrichment_group'
    x[, group := groupnames[i]]
    x[, muttype := mutnames[i]]
    gmeta[x, , on=.(group)]
}))

fwrite(tab, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
