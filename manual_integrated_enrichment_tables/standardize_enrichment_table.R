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
if (length(args) != 2) {
    stop("usage: standardize_enrichment_table.R enrichment_table.csv out.scv")
}

in.csv <- args[1]
out.csv <- args[2]

if (file.exists(out.csv))
    stop(paste('output file', out.csv, 'already exists, please delete it first'))

suppressMessages(library(data.table))


tab <- fread(in.csv)

# quantile analyses have columns QUANTILES, MIN_COVERAGE and BINSIZE (in caps)
is.bed.regions <- !('QUANTILES' %in% colnames(tab))
if (is.bed.regions) {
    tab[, QUANTILES := NA]
    tab[, MIN_COVERAGE := NA]
    tab[, BINSIZE := NA]
}

rename.column <- function(old, new) setnames(tab, old=old, new=new)

datasource <- tab[,datasource][1]
dataclass <- ''
if ('dataclass' %in% colnames(tab))
    dataclass <- tab[,dataclass][1]

if (is.bed.regions & datasource == 'roadmap' & dataclass == 'chromhmm') {
    tab <- tab[eid == "E073"]
    tab[, model := NULL]
    tab[, eid := NULL]
} else if (datasource == 'ucsc_track') {
    rename.column('track', 'dataclass')
} else if (datasource == 'gtex' & 'enrichment_group' %in% colnames(tab) & 'tissue' %in% colnames(tab)) {
    rename.column('enrichment_group', 'dataclass')
    rename.column('tissue', 'lineclass')
} else if (datasource == 'repliseq') {
    rename.column('i.celltype', 'dataclass')
    tab[, geoid := NULL]
} else if (datasource == 'encode' & dataclass == 'histone_marks') {
    tab <- tab[eid == 'E073']
    tab[, eid := NULL]
    tab[, signal_type := NULL]
    rename.column('mark', 'lineclass')
} else if (datasource == 'scatacseq') {
    tab <- tab[libid == 'librarymerged' & sample == 'merged']
    rename.column('libid', 'dataclass')
    tab[, sample := NULL]
    rename.column('i.celltype', 'lineclass')
} else if (datasource == 'scrnaseq') {
    tab <- tab[donor == 'combined' & selection == 'combined']
    tab[, donor := NULL]
    tab[, selection := NULL]
    rename.column('i.celltype', 'lineclass')
} else if (datasource == 'fragile_sites') {
} else if (datasource == 'gene_length') {
    # if we ever use a different tissue for expression, will need to represent that here
    rename.column('transcript_type', 'lineclass')
}


if (!('dataclass' %in% colnames(tab)))
    tab[, lineclass := NA] # ???: is this correct?
if (!('lineclass' %in% colnames(tab)))
    tab[, lineclass := NA]

colorder <- c(
    'group',
    'amp',
    'phenotype',
    'celltype',
    'delsig',
    'muttype',
    'QUANTILES',
    'MIN_COVERAGE',
    'BINSIZE',
    'datasource',
    'dataclass',
    'lineclass',
    'quantile',
    'pval',
    'padj',
    'fdr',
    'enr',
    'obs',
    'perm.mean',
    'perm.lb.95',
    'perm.q1',
    'perm.med',
    'perm.q3',
    'perm.ub.95',
    'enr.perm.lb.95',
    'enr.perm.q1',
    'enr.perm.med',
    'enr.perm.q3',
    'enr.perm.ub.95',
    'n.bootstraps',
    'boot.0.95.lb',
    'boot.0.95.ub',
    'enr.boot.0.95.lb',
    'enr.boot.0.95.ub')

# setcolorder() rearranges such that `colorder' named columns are at the front
# it does not drop additional unnamed columns
setcolorder(tab, colorder)
tab <- tab[, 1:length(colorder)]

fwrite(tab, file=out.csv)

if ('snakemake' %in% ls()) {
    sink()
}
