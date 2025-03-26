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
        snakemake@input['csv'],
        snakemake@params['tags'],
        snakemake@input['expomats']       # variable length list
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("usage: annotate_burdens_with_exposures.R out.csv burden_long_table.csv tag_string exposures1.csv [ exposures2.csv ... exposuresN.csv ]")
}

out.csv <- args[1]
in.csv <- args[2]
tag.string <- args[3]
exposure.csvs <- args[-(1:3)]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))

tags <- fread(gsub(pattern=';', replacement='\n', paste0(tag.string, '\n')),
    #text=strsplit(tag.string, ';')[[1]],
    header=FALSE, col.names=c('burden.source', 'burden.type', 'mut.type'))
if (nrow(tags) != length(exposure.csvs))
    stop(paste("Tag string", tag.string, "has more ;-separated fields than the number of supplied exposure CSVs", length(exposure.csvs)))

tags$file <- exposure.csvs

burden.tab <- fread(in.csv)

# burden calculations occur prior to recurrence filtering
burden.tab[, scaling.factor := burden/num.muts.unfiltered]

exposures <- rbindlist(lapply(1:nrow(tags), function(i) {
    e <- fread(tags[i,file])
    # Example file format:
    #   Samples ID8 ID9 ID11    ID21    ID83A
    #   11901106_Gy_C1  0   26  0   40  0
    #   11901106_Gy_C2  0   0   17  0   30
    #   1278-Oligo-5    1   0   0   0   0
    ret <- data.table::melt(data.table(sample=e$Samples,
        burden.source=tags[i, burden.source],
        burden.type=tags[i, burden.type],
        mut.type=tags[i, mut.type],
        num.muts.unfiltered=NA, # Signature analysis is not performed on unfiltered calls
        e[,-1]),
        measure.vars=colnames(e)[-1],  # all but samples column
        variable.name='burden.subtype',
        value.name='num.muts')
    setcolorder(ret, c('sample', 'burden.source', 'burden.type', 'burden.subtype', 'mut.type', 'num.muts.unfiltered', 'num.muts'))
    ret
}))

exposures[, burden := num.muts * burden.tab[exposures, scaling.factor, on=.(sample, mut.type)]]

fwrite(exposures, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
