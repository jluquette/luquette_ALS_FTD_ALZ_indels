library(scan2)
library(data.table)
library(ggplot2)
library(ggrepel)

meta <- fread('../metadata/sample_metadata.csv')

muts <- fread('../tables/all___FILTERED_mut___any.csv')[pass == TRUE] # do not use rescued mutations

summary.objects <- list.files(pattern='.*.rda', path='../scan2/summary_objects', full.name=T)

reload=TRUE
if (reload) {
    mapds <- rbindlist(lapply(summary.objects, function(fn) {
        print(fn)
        load(fn)
        data.table(sample=name(smry), mapd(smry))
    }))
    colnames(mapds)[3] <- 'mapd'
}

pink.peaks <- c(
     "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",  # pink peaks
     "2:Del:M:1", "3:Del:M:1", "3:Del:M:2"  # microhomology purple peaks
)

counts <- muts[,.(snvs=sum(muttype=='snv'), indels=sum(muttype=='indel'),
    pink=sum(muttype=='indel' & mutsig %in% pink.peaks),
    not.pink=sum(muttype=='indel' & !(mutsig %in% pink.peaks))), by=.(sample)]

counts <- meta[mapds[binsize==64000][counts, , on=.(sample)], , on=.(sample)]
counts <- melt(counts, measure.vars=c('snvs', 'indels', 'pink', 'not.pink'), variable.name='count.type', value.name='count')
counts <- counts[celltype == 'neuron'] # have checked with oligos: none are outliers

outliers <- unique(counts[((amp=='MDA' & mapd>=2) | (amp=='PTA' & mapd>=0.4)) & ((count.type == 'not.pink' & count >= 150) | (count.type == 'pink' & count >= 400))]$sample)

fwrite(data.table(sample=outliers), file='outliers.csv')


# 3 plots

# 1. plot with outliers included without labels, which makes the outliers easier to see
p <- ggplot(counts[count.type != 'indels'], aes(x=mapd, y=count, shape=phenotype, col=delsig)) +
    geom_point() + 
    facet_grid(count.type~amp, scale='free') +
    labs(title='Mutation counts by quality (MAPD using 64kb bins)',
        subtitle='No labels')
print(p)
dev.print(dev=pdf, file='with_outliers_no_labels.pdf')

# 2. plot with outliers included with labels to show which cells are flagged and where they
#    appear in the different mutation type plots
p <- ggplot(counts[count.type != 'indels'], aes(x=mapd, y=count, shape=phenotype, col=delsig)) +
    geom_point() + 
    geom_text_repel(data=counts[count.type != 'indels' & sample %in% outliers], aes(label=sample), min.segment.length=0, box.padding=1, max.overlaps=20) +
    facet_grid(count.type~amp, scale='free') +
    labs(title='Mutation counts by quality (MAPD using 64kb bins)',
        subtitle='Outlier cells labelled in all plots')
print(p)
dev.print(dev=pdf, file='with_outliers_with_labels.pdf')

# 3. plot with no outliers to better show distribution
p <- ggplot(counts[!(sample %in% outliers) & count.type != 'indels'],
        aes(x=mapd, y=count, shape=phenotype, col=delsig)) +
    geom_point() + 
    facet_grid(count.type~amp, scale='free') +
    labs(title='Mutation counts by quality (MAPD using 64kb bins)',
        subtitle='Outliers removed to better show distribution')
print(p)
dev.print(dev=pdf, file='no_outliers.pdf')


# 4. plot with no outliers and linear model lines showing no count vs. MAPD trend for PTA data
p <- ggplot(counts[!(sample %in% outliers) & count.type != 'indels'],
        aes(x=mapd, y=count, shape=phenotype, col=delsig, group=delsig)) +
    geom_point() + 
    geom_smooth(method='lm', se=F) +
    facet_grid(count.type~amp, scale='free') +
    labs(title='Mutation counts by quality (MAPD using 64kb bins)',
        subtitle='Outliers removed and linear model of count vs. MAPD shown')
print(p)
dev.print(dev=pdf, file='no_outliers_with_linear_models.pdf')
