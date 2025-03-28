#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['neuron'],
        snakemake@input['oligo'],
        snakemake@input['leesix'],
        snakemake@input['machado'],
        snakemake@output['csv'],
        snakemake@output['pdf'],
        snakemake@output['svg']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop('usage: fig2_panel_a.R neuron_spectrum.csv oligo_spectrum.csv leesix_spectrum.csv machado_sbs.csv out.csv out.pdf out.svg')
}


neuron.csv <- args[1]
oligo.csv <- args[2]
leesix.csv <- args[3]
machado.csv <- args[4]
out.csv <- args[5]
out.pdf <- args[6]
out.svg <- args[7]

for (f in c(out.csv, out.pdf, out.svg)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(svglite))


neuron <- fread(neuron.csv)
neuron.spectrum <- neuron$Spectrum
oligo <- fread(oligo.csv)
oligo.spectrum <- oligo$Spectrum
machado <- fread(machado.csv)
machado <- machado$SBSblood / sum(machado$SBSblood)
leesix <- fread(leesix.csv)
leesix <- leesix$HSPC_spectrum / sum(leesix$HSPC_spectrum)

cossim <- function(a, b) sum(a*b) /(sqrt(sum(a^2))*sqrt(sum(b^2)))

o.cossim <- sapply(list(neuron.spectrum, leesix, machado), function(x) cossim(oligo.spectrum, x))
n.cossim <- sapply(list(oligo.spectrum, leesix, machado), function(x) cossim(neuron.spectrum, x))

d <- data.table(MutType=neuron$MutType, Neurons=neuron, Oligo=oligo,
    LeeSix2018=leesix, Machado2020=machado)
fwrite(d, file=out.csv)

devs=list(pdf, svglite)
outs=c(out.pdf, out.svg)
for (i in 1:2) {
    devs[[i]](width=3, height=4, pointsize=5, file=outs[i])
    layout(1:4)
    par(mar=c(1,4,2,1))
    plot.sbs96(x=0, spectrum=neuron.spectrum, main='Neurons')
    legend('topright', legend=paste(c('Oligo', 'Lee-Six HSPC', 'Machado blood'), round(n.cossim,3)),
        bty='n', bg='white', title='Cosine sim.')
    plot.sbs96(x=0, spectrum=oligo.spectrum, main='Oligodendrocytes')
    legend('topright', legend=paste(c('Neuron', 'Lee-Six HSPC', 'Machado blood'), round(o.cossim,3)),
        bty='n', bg='white', title='Cosine sim.')
                                      
    plot.sbs96(x=0, spectrum=leesix, main='HSPC spectrum (Lee-Six et al. 2018)')
    plot.sbs96(x=0, spectrum=machado, main='SBSblood (Machado et al. 2022)')

    dev.off()
}

if ('snakemake' %in% ls()) {
    sink()
}
