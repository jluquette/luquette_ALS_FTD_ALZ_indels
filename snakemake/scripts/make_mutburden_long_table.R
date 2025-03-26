#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['csv'],
        snakemake@input['meta'],
        snakemake@input['filtered_muts'],
        snakemake@input['objects']       # variable length list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("usage: make_mutburden_long_table.R out.csv metadata.csv filtered_mutations.csv scan2_object1.rda [ scan2_object2.rda .. scan2_objectN.rda]")
}

out.csv <- args[1]
meta <- args[2]
filtered.muts.csv <- args[3]
object.rdas <- args[-(1:3)]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(scan2))

meta <- fread(meta)
filt.muts <- fread(filtered.muts.csv)

burdens <- rbindlist(lapply(object.rdas, function(object.path) {
    o <- load.summary(paths=object.path)
    data.table(sample=name(o),
        burden.source='scan2',
        burden.type='extrapolation',
        burden.subtype='uncorrected',
        mut.type=c('snv', 'indel'),
        num.muts.unfiltered=c(nrow(passing(o, muttype='snv')), nrow(passing(o, muttype='indel'))),
        # Filtering refers to recurrence filters that are applied after SCAN2's calls.
        # These filters are not recorded in the objects.
        num.muts=filt.muts[sample == name(o) & rescue == FALSE, c(sum(muttype=='snv'), sum(muttype=='indel'))],
        burden=c(mutburden(o, muttype='snv'), mutburden(o, muttype='indel')))
}))

meta <- meta[burdens, on=.(sample)]

fwrite(meta, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
