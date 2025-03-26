#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@output['svg'],
        snakemake@output['pdf'],
        snakemake@output['burdens'],
        snakemake@output['models'],
        paste(snakemake@params[['colors']], collapse=','),
        paste(unlist(snakemake@params['group_names']), snakemake@input, sep='=')
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    cat('`colors` is a comma separated string of R recognized color strings IN THE SAME ORDER AS THE burden_files.\n')
    cat("`colors=table` to use the colors already present in the burden files.\n")
    stop("usage: model_aging_rates.r out.svg out.pdf out.combined_burdens.csv out.models.csv colors group1_name=burden_file1.csv [ group2_name=burden_file2.csv ... groupN_name=burden_fileN.csv ]")
}

outsvg <- args[1]
outpdf <- args[2]
outcsv <- args[3]
outmodelcsv <- args[4]
color.string <- args[5]
burden.csvs <- args[-(1:5)]

using.table.colors=TRUE
if (color.string != "table") {
    colors <- unlist(strsplit(color.string, split=',')[[1]])
    using.table.colors <- FALSE
    print(colors)
} else {
    colors <- color.string
    print("colors=table, not overriding table colors")
}

if (file.exists(outsvg))
    stop(paste('output file', outsvg, 'already exists, please delete it first'))
if (file.exists(outpdf))
    stop(paste('output file', outpdf, 'already exists, please delete it first'))
if (file.exists(outcsv))
    stop(paste('output file', outcsv, 'already exists, please delete it first'))
if (file.exists(outmodelcsv))
    stop(paste('output file', outmodelcsv, 'already exists, please delete it first'))

group.names <- sapply(strsplit(burden.csvs, '='), head, 1)
burden.csvs <- sapply(strsplit(burden.csvs, '='), function(x) paste(x[-1], collapse='='))

suppressMessages(library(data.table))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(svglite))


burdens <- setNames(lapply(1:length(burden.csvs), function(i) {
    in.csv <- burden.csvs[i]
    burden <- fread(in.csv)
    burden[, group := group.names[i]]
    if (!using.table.colors) {
        cat(paste0('group=', group.names[i], ': overiding table color=', burden$color[1], ' with color=', colors[i], '\n'))
        burden[, color := colors[i]]
    }
    if (!('uncorrected.genome.burden' %in% colnames(burden))) {
        cat(paste0('group=', group.names[i], ': adding copyover uncorrected.genome.burden column to table\n'))
        burden[, uncorrected.genome.burden := genome.burden]
    }
    burden
}), group.names)
print(burdens)

# harmonize column names/order before rbindlist
# unique maintains order
all.col.names <- unique(unlist(lapply(burdens, colnames)))
burdens <- lapply(burdens, function(b) setcolorder(b, all.col.names))

combined.burdens <- rbindlist(burdens)
print(combined.burdens)

burdens <- c(burdens, `All groups`=list(combined.burdens))
print(burdens)
models <- lapply(burdens, function(dt) {
print(table(dt$group))
print(dt)
    if (length(unique(dt$group)) == 1) {
        if (length(unique(dt$donor)) == 1)
            model <- lm(genome.burden ~ age, data=dt) #[outlier == "NORMAL"])
        else
            model <- lmer(genome.burden ~ age + (1|donor), data=dt) #[outlier == "NORMAL"])
    } else {
        # for the combined model
        if (length(unique(dt$donor)) == 1)
            model <- lm(genome.burden ~ age*group, data=dt) #[outlier == "NORMAL"])
        else
            model <- lmer(genome.burden ~ age*group + (1|donor), data=dt) #[outlier == "NORMAL"])
    }
    ci <- confint(model)
    colnames(ci) <- paste0('95% CI ', c('lower', 'upper'))
    list(model=model, ci=ci)
})
print(models)

figwidth=4*3   # add two extra panels (one left, one right):
               # left: label to denote with or without outliers; right legend
figheight=4*2  # two rows: top row with outliers included as Xs, bottom no outliers
devs <- list(svglite, pdf)
outs <- c(outsvg, outpdf)
for (i in 1:2) {
    devs[[i]](file=outs[i], width=figwidth, height=figheight)

    layout(matrix(1:6,nrow=2,byrow=T))

    #for (include.outliers in c(F,T)) {
    for (include.outliers in c(T)) {
        plot(1, pch=NA, bty='n', xlab='', ylab='', xaxt='n', yaxt='n')
        legend('center', bty='n', legend=paste("Outliers included:", include.outliers))

        # any donor can be used
        # [-length(models)] - don't use the combined model
        # XXX: no longer works due to mixing lmer() and lm()
        #maxpred <- sapply(models[-length(models)], function(m)
            #predict(m$model, data.frame(age=max(combined.burdens$age), donor=m$model@frame$donor[1]))
        #)
        zzz <- combined.burdens
        if (!include.outliers)
            zzz <- zzz[outlier == 'NORMAL']

        # make ylim big enough to contain all data points and lines
        #ylim <- c(0, max(zzz$genome.burden, maxpred, na.rm=TRUE))
        ylim <- c(0, max(zzz$genome.burden, na.rm=TRUE))
        plot(zzz[, .(age, genome.burden)],
            col=zzz$color, pch=ifelse(zzz$outlier == 'NORMAL', 17, 4),
            bty='l', cex=1,
            ylim=ylim, ylab='Autosomal mutation burden', xlab='Age')
        abline(h=0, col='grey')
        # These curve() calls don't plot below 0.
        for (i in 1:(length(models)-1)) {
            m <- models[[i]]
            curve(coef(summary(m$model))[1] + coef(summary(m$model))[2]*x,
                # don't use zzz here instead of combined.burdens: keep age axis identical
                from=0, to=2*max(combined.burdens$age), lwd=2,
                col=burdens[[i]]$color[1], add=T)
        }
    
        # make new plot for legend
        plot(1, pch=NA, bty='n', xlab='', ylab='', xaxt='n', yaxt='n')
        legend('center', lwd=2, pch=17, bty='n',
            col=sapply(burdens[-length(burdens)], function(b) b$color[1]),
            legend=names(burdens)[-length(burdens)])
    }
}

fwrite(combined.burdens, file=outcsv)

for (i in 1:length(models)) {
    m <- models[[i]]
    print(names(models)[i])
    print(summary(m$model))
}

all.models <- rbindlist(lapply(1:length(models), function(i) {
    m <- models[[i]]
    coefs <- coef(summary(m$model))
    data.table(Group=names(models)[i],
        Formula=deparse(formula(m$model)),
        Variable=rownames(coefs),
        cbind(coefs, m$ci[rownames(coefs),]))
    }),
    # XXX: after mixing lm() and lmer(), there are different numbers of columns
    # no longer know if they can be rbind()ed, so this table might not be valid
    # anymore.
    fill=TRUE)
fwrite(all.models, file=outmodelcsv)

if ('snakemake' %in% ls()) {
    sink()
}
