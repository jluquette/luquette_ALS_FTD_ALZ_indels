#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['exps'],
        snakemake@input['meta'],
        snakemake@output['csv'],
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("usage: residualize_burden.NEW.R out.csv sample_metadata.csv combined_burden_long_table.no_metadata.csv")
}

outcsv <- args[1]
inmeta <- args[2]
burdens.csv <- args[3]

for (f in c(outcsv)) {
    if (file.exists(f)) {
        stop(paste('output file', f, 'already exists, please delete it first'))
    }
}

suppressMessages(library(data.table))

meta <- fread(inmeta)
burdens <- fread(burdens.csv)
burdens <- meta[burdens,,on=.(sample)]  # annotate burdens with age, phenotype, etc

# For now, set burden.subtype=mut.type for SCAN2 total burden estimates
burdens[burden.source=='scan2', burden.subtype := mut.type]

# Derive linear model from control data to show how non-control data differs
# A different linear model is needed for each amplification type (i.e., PTA and MDA)
# and burden type (e.g., scan2 total burden, denovo sig, COSMIC sigs) and
# burden subtype (e.g., snvs, indels, SBSA, SBS1, IDA, ID1, ...)
cat("Residualizing burdens against: phenotype=='Control' & celltype=='neuron'\n")
controls <- burdens[phenotype == 'Control' & celltype == 'neuron']
for (bsrc in unique(burdens$burden.source)) {
for (btype in unique(burdens$burden.type)) {
for (bsubtype in unique(burdens$burden.subtype)) {
    for (this.amp in c('MDA', 'PTA')) {
        this.data <- controls[amp == this.amp & burden.source == bsrc & burden.type == btype & burden.subtype == bsubtype,]
        # not all combinations of unique (source, type, subtype, amp) exist in the data.
        # nrow(this.data)=0 when this is the case
        if (nrow(this.data) == 0)
            next
print(this.data)
cat('------------------------------------ FIX ME ----------------------------------------\n')
cat(' Currently using non-signature B corrected burdens for MDA cells!!! \n')
cat('------------------------------------ FIX ME ----------------------------------------\n')
        cat(bsrc, btype, bsubtype, this.amp, '------------------------------------------------------- \n')
        # Compute which entries are used once
        # IMPORTATNT!! `inclusion` applies to the burdens table, not controls
        inclusion <- burdens[, amp == this.amp & burden.source == bsrc & burden.type == btype & burden.subtype == bsubtype]
#print(this.data)
        # First model uses all points, including extreme outliers
        m1 <- lm(burden ~ age, data=this.data)

        # Determine outliers by applying Tukey method to residuals
        r1 <- resid(m1)
        cat("Ignoring", nrow(this.data)-length(r1), " NAs / ", nrow(this.data), " total observations\n")
        q1 <- quantile(r1, 0.25)
        q3 <- quantile(r1, 0.75)
        ub <- q3 + (q3-q1)*1.5
        lb <- q1 - (q3-q1)*1.5
        cat('tukey on residuals: lower', lb, 'upper', ub, '\n')
        this.data[!is.na(burden), tukey.outlier := r1 < lb | r1 > ub]
print(this.data)
        # Record for later in case it's useful
        burdens[phenotype == 'Control' & celltype == 'neuron' & inclusion &
            #amp == this.amp & burden.source == bsrc & burden.type == btype & burden.subtype == bsubtype &
            #& burden.type == mt & amp == this.amp &
            !is.na(burden), tukey.outlier := r1 < lb | r1 > ub]

        # Model only on this control data, excluding outliers
        m <- lm(burden ~ age, data=this.data[tukey.outlier == FALSE])
        r2 <- resid(m)
        resid.mean <- mean(r2)
        resid.sd <- sd(r2)
        cat("Final model (outliers removed)\n")
        print(summary(m))
        print(summary(r2))
        cat("resid mean:", resid.mean, '\n')
        cat("resid sd:", resid.sd, '\n')

        # Predict on all data
        burdens[inclusion, expected.burden := predict(m, newdata=burdens[inclusion])]
        #burdens[burden.type == mt & amp == this.amp, resid.burden := burden - expected.burden]
        burdens[inclusion, resid.burden := burden - expected.burden]
        #burdens[burden.type == mt & amp == this.amp, std.resid.burden := (resid.burden - resid.mean) / resid.sd]
        burdens[inclusion, std.resid.burden := (resid.burden - resid.mean) / resid.sd]
    }
}
}
}

cat("Writing long table to", outcsv, "...\n")
fwrite(burdens, file=outcsv)

if ('snakemake' %in% ls()) {
    sink()
}
