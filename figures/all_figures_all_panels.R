# STEPS FOR PRODUCING THE FINAL COLLECTED TABLES OF BURDENS
#
# 1. Create (channel X sample) matrices for SigProfilerExtractor
#   -> mutsigs/matrices/all/matrix.{SBS96,ID83,ID83_corrected}.txt
#   1a. Run SigProfilerMatrixGenerator (n=544 cells out of n=547. Does not include MDA FTD, e.g.)
#   1b. Correct ID83 matrices for SCAN2 channel-specific sensitivity
# 2. Run SigProfilerExtractor on each matrix
#   -> mutsigs/sigprofilerextractor/all/{SBS96,ID83,ID83_corrected}/Most_Stab_Sigs/{signatures,exposures}.csv
#   1a. Parallelized and code to extract the most stable signature (can be different from
#       SigProfilerExtractor's selected solution). Selected signature is copied into Most_Stab_Sigs
# 3. Create final long table of signature burdens (sample, signature types) [not a matrix]
#   -> sigburden_long_table.csv
#   3a. code_from_laptop/small_for_download/make_sigburden_longtable.sh
#       . snakemake/scripts/annotate_burdens_with_exposures.R combines SNV and indel
#         scaling factors with SBS and ID signature counts.
# 4. Combine signature burdens with SNV and indel burdens
#   -> combined_burdens_denovo_cosmic.no_metadata.csv
#   4a. Remove metadata columns from aging_rates/mutation_burdens_long_table.csv.
#       . This is just a `cut` to remove the 10 metadata columns (excluding sample).
# 5. Fit control models to compute expected burdens
#   5a. residualize_burden.NEW.R
#       . This computes a standard linear model (lm, not lmer) on a subset of control
#         PTA neurons. The subset is created by running a preliminary lm() and excluding
#         neurons identified by Tukey's outlier procedure.
#       . Expected values for all PTA cells (of the same type: neuron or oligo) are then
#         computed from the linear model without regard to phenotype.

########################################################################################
#
# Code below copied from code_from_laptop/small_for_download/zinan_cshl_talk/code.R
#
########################################################################################

library(ggplot2)
library(forcats)
library(ggsignif)
library(data.table)
library(ggpmisc)  # for p-values on correlation plots
library(scan2)    # for signature plotting
library(patchwork)

base.text=10
large.text=12
signif.size=1.6 # size of the p-value labels on significance tests. not sure what the unit is

base.theme <- list(
    theme_classic(),
    theme(axis.title=element_text(size=large.text),
          axis.text=element_text(size=base.text),
          axis.line=element_line(linewidth=0.25),
          axis.ticks=element_line(linewidth=0.25),
          strip.text=element_text(size=large.text),
          strip.background=element_blank(),
          legend.title=element_text(size=large.text),
          legend.text=element_text(size=base.text)) #,
          #aspect.ratio=1)
)

#phenotype.colors <- c(Control='#cccccc', ALS='#D66929', FTD='#F3BB71', `ALS/FTD`='#F3BB71', Alzheimers='#B1ABCF') #, na.value='#aaaaaa')

# AD and Alzheimers refer to the same thing. Eventually standardize on AD.
phenotype.colors <- c(Control='#cccccc', ALS='#DF6236', FTD='#F4C045', `ALS/FTD`='#F4C045', Alzheimers='#6C559F', AD='#6C559F') #, na.value='#aaaaaa')


ggtheme <- list(
    base.theme,
    scale_colour_manual(name='phenotype', values=phenotype.colors), #, na.value='grey'),
    scale_fill_manual(name='phenotype', values=phenotype.colors) #, na.value='grey')
)


meta <- fread('metadata/sample_metadata.csv')
meta[, phenotype := factor(phenotype, levels=c('Control', 'ALS', 'FTD', 'Alzheimers'))]
donor.order <- meta[, .(donor, age), by=donor][order(age)][!duplicated(donor)]$donor
meta[, donor := factor(donor, levels=donor.order)]
sample.order <- meta[, .(sample, age), by=sample][order(age, sample)][!duplicated(sample)]$sample
meta[, sample := factor(sample, levels=sample.order)]


groups <- fread('metadata/groups.csv')
no.mapd.indel.outliers <- groups[grep(x=group, pattern='no_mapd_indel_outliers'), unique(sample)]

########################################################################################
#
# END COPIED CODE
#
########################################################################################


burdens <- fread('aging_rates/mutation_burdens_long_table.RESIDUALIZED.csv')

# recode snv,indel -> ordered(SNV, Indel)
burdens[, mut.type := forcats::fct_shift(forcats::fct_recode(mut.type, SNV='snv', Indel='indel'))]
# copy recodings done to meta above
burdens[, phenotype := factor(phenotype, levels=c('Control', 'ALS', 'FTD', 'Alzheimers'))]
burdens[, donor := factor(donor, levels=donor.order)]
burdens[, sample := factor(sample, levels=sample.order)]

#
# 1. Show BA6 and PFC neurons have similar SNV and indel burdens, after correcting for age.
#
# Without age limit (to ensure comparable data), p-value becomes marginally significant (~0.02)
# This could mean the data doesn't reach significance because of N or that aging correction
# is not fully effective. The latter is certainly true to some extent.
ggplot(burdens[age >= 50 & sample %in% no.mapd.indel.outliers & phenotype=='Control' & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=tissue.origin, y=burden/expected.burden, fill=phenotype)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.2, size=1/2) +
    geom_signif(comparisons=list(c('BA6', 'PFC')), margin_top=-0.1, textsize=3.25) +
    geom_hline(yintercept=1, linetype=2) +
    facet_grid( mut.type ~ .) + #, scale='free_y') +
    xlab('Brain region') + ylab("Mutation burden (obs / exp)") +
    ggtheme + theme(aspect.ratio=2) + expand_limits(y=c(0, 2)) +
    labs(title='Control neuron burdens', subtitle='PTA, no MAPD/indel outliers, age>50')
dev.print(dev=pdf, file='ba6_vs_pfc_control_neuron_burden_comparison.age_gt_50.pdf')



#
# 2. Per individual burden boxes - show that in some individuals, intra-individual burden can vary greatly
#
# Need to split by BA9 vs. BA6?
ggplot(burdens[sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=donor, y=burden, fill=tissue.origin)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(size=1/5, width=0.2) +
    facet_grid(forcats::fct_shift(forcats::fct_recode(mut.type, SNV='snv', Indel='indel')) ~ phenotype, scale='free', space='free_x') +
    facet_grid(mut.type ~ phenotype, scale='free', space='free_x') +
    xlab('Individual ID') + ylab("Mutation burden (log-scale)") +
    base.theme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + scale_y_log10()
    #+ scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2))
dev.print(dev=pdf, file='burdens_per_individual.pdf')



#
# 3. Aging trend lines and points - trends are only shown for controls for comparison
#
ctrl.burdens <- burdens[phenotype == 'Control' & amp == 'PTA' & celltype == 'neuron']
ctrl.burdens[, phenotype := NULL]
ggplot(burdens[phenotype != 'Control' & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=age, y=burden, col=phenotype)) +
    geom_point(data=ctrl.burdens, aes(age, burden, col='Control'), size=1/2) +
    geom_smooth(method='lm', se=F, data=ctrl.burdens, aes(age, burden, col='Control'), size=1/2) +
    geom_point(size=1.5) +
    facet_grid(mut.type~phenotype, scales='free_y') +
    xlab('Age') + ylab('Somatic mutation burden') + ggtheme + theme(aspect.ratio=1) +
    labs(title='Mutation burdens with control trend lines', subtitle='MAPD/indel outliers removed')
dev.print(dev=pdf, file='mutation_burdens_with_control_trend_lines.NO_MAPD_INDEL_OUTLIERS.pdf')


#
# 4. Boxplots comparing mutation residual-adjusted age burdens between phenotypes
#
ggplot(burdens[age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=phenotype, y=burden/expected.burden, fill=phenotype)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.25, size=1/2) +
    geom_hline(yintercept=1, linetype=2, linewidth=0.25) +
    geom_signif(comparisons=list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')), step_increase=1/12, margin_top=-0.2) +
    scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2)) +
    xlab('') + ylab('Mutation burden (obs / exp, log-scale)') +
    labs(title='Burden residual', subtitle='MAPD/indel outliers removed, age>50 required') +
    expand_limits(y=c(0.5, 2.25)) +
    facet_grid(mut.type~., scale='free_y') +
    ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2))
dev.print(dev=pdf, file='residualized_mutation_burden_ratios.by_phenotype.NO_MAPD_INDEL_OUTLIERS_AGE_GT_50.pdf')



#
# 5. Boxplots comparing mutation residual-adjusted age burdens between TDP-43 statuses
#
ggplot(burdens[age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA' & phenotype %in% c('ALS', 'FTD')], aes(x=fct_relevel(forcats::fct_cross(phenotype, fct_shift(forcats::fct_recode(tdp43, `TDP43-`='Depleted', `TDP43+`='Normal')), sep=' '), c('ALS TDP43+', 'ALS TDP43-', 'FTD TDP43+', 'FTD TDP43-')), y=burden/expected.burden, fill=phenotype)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.25, size=1/2) +
    geom_hline(yintercept=1, linetype=2, linewidth=0.25) +
    geom_signif(comparisons=list(c('ALS TDP43-', 'ALS TDP43+'), c('FTD TDP43-', 'FTD TDP43+')), textsize=2.5, size=1/4, margin_top=0.00, y_position=log10(25)) + #, position=position_nudge(y=-0.2)) +
    xlab('') + ylab('Mutation burden (obs / exp, log-scale)') +
    scale_y_log10() + expand_limits(y=c(0.75, 30)) +
    labs(title='Burden residual', subtitle='MAPD/indel outliers removed, age>50 required') +
    facet_grid(~mut.type) + #, scale='free_y', ncol=2) +
    ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2), aspect.ratio=3/2)
dev.print(dev=pdf, file='residualized_mutation_burden_ratios.by_tdp43.NO_MAPD_INDEL_OUTLIERS_AGE_GT_50.REORDERED.pdf')
    #geom_signif(comparisons=list(c('ALS TDP43-', 'ALS TDP43+'), c('FTD TDP43-', 'FTD TDP43+')), margin_top=-0.1, textsize=2.5, size=1/4) +
    #scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2)) + #, expand=expansion(add=1.10)) +



#
# 6a. Boxplots comparing denovo SBS96 burdens (residual-adjusted age)
#
resids <- meta[fread('combined_residualized_burdens.csv'), , on=.(sample)]
resids[, obs.exp := (1+burden)/(1+expected.burden)]
ggplot(resids[burden.source=='denovo' & burden.type=='SBS96' & age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=phenotype, y=(1+burden)/(1+expected.burden), fill=phenotype)) +
    geom_boxplot(outlier.shape=NA, col=1, linewidth=0.3) +
    geom_jitter(col=1, size=1/2) +
    facet_wrap(~burden.subtype, ncol=3) +
    ggtheme +
    xlab('') + ylab('Signature exposure (obs / exp)') +
    guides(fill=guide_legend(title='Phenotype')) +
    geom_hline(yintercept=1, linewidth=0.5) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), aspect.ratio=2.0) +
    geom_signif(comparisons=list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')),
        step_increase=1/10, margin_top=-0.1, textsize=2.75, size=0.15) +
    labs(title='De novo SBS96 burdens by phenotype', subtitle='MAPD/indel outliers removed, age>50')
    #scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2))
dev.print(dev=pdf, 'denovo_sbs96_residualized_burdens_age_matched.NO_MAPD_OUTLIERS_AGE_GT_50.pdf')

#
# 6b. Age vs. denovo SBS96 burdens
#
ctrl.resids <- resids[phenotype=='Control' & amp == 'PTA' & celltype=='neuron']
ctrl.resids[, phenotype := NULL]
p1 <- ggplot(resids[phenotype != 'Control' & burden.source=='denovo' & burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[1:3]) & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=age, y=burden, col=phenotype)) +
    geom_point(data=ctrl.resids[burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[1:3])], aes(age, burden, col='Control'), size=1/2) +
    geom_smooth(method='lm', se=F, data=ctrl.resids[burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[1:3])], aes(age, burden, col='Control'), size=1/2) +
    geom_point(size=1) +
    facet_grid(burden.subtype~phenotype, scales='free_y') +
    xlab('Age') + ylab('De novo signature burden') + ggtheme + theme(aspect.ratio=1)
p2 <- ggplot(resids[phenotype != 'Control' & burden.source=='denovo' & burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[4:6]) & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=age, y=burden, col=phenotype)) +
    geom_point(data=ctrl.resids[burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[4:6])], aes(age, burden, col='Control'), size=1/2) +
    geom_smooth(method='lm', se=F, data=ctrl.resids[burden.type=='SBS96' & burden.subtype %in% paste0('SBS96', LETTERS[4:6])], aes(age, burden, col='Control'), size=1/2) +
    geom_point(size=1) +
    facet_grid(burden.subtype~phenotype, scales='free_y') +
    xlab('Age') + ylab('De novo signature burden') + ggtheme + theme(aspect.ratio=1)
(p1|p2)+plot_layout(guides='collect')
dev.print(dev=pdf, file='denovo_sbs96_burdens_vs_age_with_control_trend_lines.NO_MAPD_INDEL_OUTLIERS.pdf')



#
# 7a. Boxplots comparing denovo ID83 burdens (residual-adjusted age)
#
ggplot(resids[burden.source=='denovo' & burden.type=='ID83_corrected' & age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=phenotype, y=pmax(1/4, (1+burden)/(1+expected.burden)), fill=phenotype)) +
    geom_boxplot(outlier.shape=NA, col=1, linewidth=0.3) +
    geom_jitter(col=1, size=1/2) +
    facet_wrap(~burden.subtype, ncol=3) +
    ggtheme +
    xlab('') + ylab('Signature exposure (obs / exp)') +
    guides(fill=guide_legend(title='Phenotype')) +
    geom_hline(yintercept=1, linewidth=0.5) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), aspect.ratio=2.0) +
    geom_signif(comparisons=list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')),
        step_increase=1/10, margin_top=-0.1, textsize=2.75, size=0.15) +
    labs(title='De novo ID83 burdens by phenotype', subtitle='MAPD/indel outliers removed, age>50, floor set=0.25') +
    scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2))
dev.print(dev=pdf, 'denovo_id83_residualized_burdens_age_matched.NO_MAPD_OUTLIERS_AGE_GT_50.pdf')


#
# 7b. Age vs. de novo ID83 burdens
#
ggplot(resids[phenotype != 'Control' & burden.source=='denovo' & burden.type=='ID83' & burden.subtype %in% c('ID83A', 'ID83B') & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=age, y=burden, col=phenotype)) +
    geom_point(data=ctrl.resids[burden.type=='ID83' & burden.subtype %in% c('ID83A', 'ID83B')], aes(age, burden, col='Control'), size=1/2) +
    geom_smooth(method='lm', se=F, data=ctrl.resids[burden.type=='ID83' & burden.subtype %in% c('ID83A', 'ID83B')], aes(age, burden, col='Control'), size=1/2) +
    geom_point(size=1) +
    facet_grid(burden.subtype~phenotype, scales='free_y') +
    xlab('Age') + ylab('De novo signature burden') + ggtheme + theme(aspect.ratio=1)
dev.print(dev=pdf, file='denovo_id83_burdens_vs_age_with_control_trend_lines.NO_MAPD_INDEL_OUTLIERS.pdf')


#
# 8. Boxplots comparing pink peak classifications by 2bp deletion amount
# NO LONGER USED IN MAIN FIGURE - but gives important plot structure for Chris'
# favorite plot: the boxplot of %2bp dels with point colors by ID-A high status
#
class.tab <- muts[age>50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & amp == 'PTA', .(donor=donor[1], sex=sex[1], age=age[1], phenotype=phenotype[1], amp=amp[1], tissue.origin=tissue.origin[1], tdp43=tdp43[1], celltype=celltype[1], delsig=fct_collapse(delsig[1], Negative='negative', Positive=c('normal', 'high')), del2=sum(grepl(pattern='2:Del', x=mutsig)), indel=sum(muttype=='indel')), by=.(sample)][indel < 20, delsig := 'Excluded']
class.tab[celltype == 'oligo', phenotype := 'Control oligodendrocyte']
class.tab[, phenotype := fct_shift(phenotype, n=-1)]
class.tab[, delsig.inclusion := fct_collapse(delsig, 'Included'=c('Negative', 'Positive'))]

# Fisher 2x2 contingency table test to see if phenotype is independent of amount of
# pink peak positive cells.
# Have to manually compute these because the y-axis values on the plot are not what
# we are testing, but those are the only things supplied to geom_signif.
fishertests <- rbindlist(lapply(list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')), function(comp) {
    x <- class.tab[delsig.inclusion=='Included' & phenotype %in% comp]
    ft <- fisher.test(x=fct_drop(x$phenotype), y=x$delsig)
    ci <- ft$conf.int
    data.frame(pheno1=comp[1], pheno2=comp[2],
        pheno1.n=sum(x$phenotype==comp[1]), pheno1.positive=sum(x[,phenotype==comp[1] & delsig=='Positive']),
        pheno2.n=sum(x$phenotype==comp[2]), pheno2.positive=sum(x[,phenotype==comp[2] & delsig=='Positive']),
        p.value=ft$p.value, odds.ratio=ft$estimate, ci.lower=ci[1], ci.upper=ci[2], ci.level=attr(ci, 'conf.level'))
}))
fishertests[, pheno1.frac := pheno1.positive/pheno1.n]
fishertests[, pheno2.frac := pheno2.positive/pheno2.n]

ggplot(class.tab, aes(x=phenotype, y=del2/indel, col=delsig, group=phenotype)) +
    geom_boxplot(outlier.shape=NA, col=1, linewidth=0.3) +
    geom_jitter(width=0.25, size=1, aes(shape=delsig.inclusion)) +
    base.theme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab('') + ylab('Fraction 2 bp deletions') +
    scale_colour_manual(name='delsig', values=c(Excluded='#000000', Negative='#000000', Positive='#CE5E70')) +
    scale_shape_manual(name='delsig.inclusion', values=c(Excluded=1, Included=16)) +
    labs(title='Disease deletion signature', subtitle='Fisher 2x2 tests of phenotype vs. #positive') +
    expand_limits(y=c(0,1)) +
    geom_signif(comparisons=list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')),
        step_increase=1/10, margin_top=-0.25, textsize=2.75, size=0.15, annotation=sprintf('%0.2g', fishertests$p.value))
dev.print(dev=pdf, file='pink_peak_class_rates_boxplot_fisher_test_pvals.NO_MAPD_OUTLIERS_AGE_GT_50.pdf')




#
# 9. ID83A vs. SBS{A,B,C...} correlation plots
# 9a. supplementary version with all 6 de novo signatures. main version is only SBS96B
#

# melted such that every row is a sample with ID83A and B values and ONE of the six SBS96_X burden values
# so there are 6 x (#samples) rows.
recast.resids <- melt(dcast(resids[burden.source=='denovo' & (burden.type=='ID83_corrected' | burden.type == 'SBS96') & age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'],  donor + sex + age + phenotype + source + sample  +  amp + tissue.origin + tdp43 + celltype +  delsig ~ burden.subtype, value.var='burden'), measure.vars=paste0('SBS96', LETTERS[1:6]), variable.name='sbs.sig')

ggplot(recast.resids[sbs.sig=='SBS96B'], aes(x=1+ID83A, y=value, col=phenotype)) +
    geom_point(size=1/2) +
    geom_smooth(method='lm', se=FALSE, linewidth=2/3) +
    ggpmisc::stat_fit_glance(method='lm', aes(label=paste('P =', round(after_stat(p.value), digits=3)))) +
    scale_x_log10() +
    ggtheme + theme(aspect.ratio=1) +
    xlab('ID83A exposure (log-scale)') + ylab('SBS96B signature exposure') +
    labs(title='Correlations between signature exposures', subtitle='No MAPD/indel outliers, no age<50 controls, P-values are lm coefficient Pr(>|t|) values, x=1+ID83A')
dev.print(dev=pdf, file='correlation_sbs96_vs_id83.pdf')

ggplot(recast.resids, aes(x=1+ID83A, y=value, col=phenotype)) +
    geom_point(size=1/2) +
    geom_smooth(method='lm', se=FALSE, linewidth=2/3) +
    ggpmisc::stat_fit_glance(method='lm', aes(label=paste('P =', round(after_stat(p.value), digits=3))), size=2.5) +
    facet_wrap(~sbs.sig, ncol=3, scale='free_y') +
    scale_x_log10() +
    ggtheme + theme(aspect.ratio=1) +
    xlab('ID83A exposure (log-scale)') + ylab('De novo SBS signature exposure') +
    labs(title='Correlations between signature exposures', subtitle='No MAPD/indel outliers, no age<50 controls, P-values are lm coefficient Pr(>|t|) values, x=1+ID83A')
dev.print(dev=pdf, file='correlation_sbs96_vs_id83.SUPPLEMENTARY.pdf') #, width=6, height=3)


#
# 9. Dinucleotide deletion freqs
#
# to collapse reverse complement-equivalent dinucs
revcomp <- c(
    CG='CG', TA='TA', AT='AT', GC='GC',

    # equiv pairs (rev comp means equivalent deletion)
    CC='CC', GG='CC',
    TT='TT', AA='TT',
    CA='CA', TG='CA',
    CT='CT', AG='CT',
    TC='TC', GA='TC',
    GT='GT', AC='GT'
)
# dinucs that we will read on the opposite strand
will.revcomp = c('GG', 'AA', 'TG', 'AG', 'GA', 'AC')


# UPDATE: 03/14/2025: add MDA data to the table, we can filter later depending on whether
# we want to show MDA or not.
# dinucs in hg19
dinucs <- fread('dinucs.txt', col.names=c('dn', 'count'))[,.(count=sum(count)),by=revcomp[dn]]
colnames(dinucs)[1] <- 'dn'
dinucs[, phenotype := 'Ref_hg19']
dinucs[, amp := NA]

dinuc.tab <- rbind(melt(dcast(muts[nchar(refnt)-nchar(altnt)==2 & muttype=='indel' & age>50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype=='neuron'], amp+phenotype ~ revcomp[dn]), id.vars=c('amp','phenotype')), dinucs[, .(amp, phenotype, variable=dn, value=count)])

norm.by.ref <- dinucs[dinuc.tab,,on=.(dn=variable)][,.(dn, obs=value, obs.frac=value/sum(value)),by=.(i.amp,i.phenotype)][,.(amp=i.amp, phenotype=i.phenotype, dn, obs, obs.frac)]


ggplot(norm.by.ref[amp != 'MDA' & phenotype != 'Ref_hg19'], aes(x=dn, y=obs.frac, fill=phenotype, group=interaction(amp, phenotype))) +
    geom_col(position='dodge', col=1) +
    ggtheme +
    xlab('Dinucleotide deletion') + ylab('Fraction of 2 bp deletions') +
    labs(title='Deleted dinucleotides', subtitle='No MAPD/indel outliers, age>50')
dev.print(dev=pdf, file='dinuc_deletion_fractions.NO_MAPD_INDEL_OUTLIERS_AGE_GT_50.pdf')



#
# 10. 10-mer logo plots
#
#devtools::install_github("omarwagih/ggseqlogo")  # if not installed
library(ggseqlogo)
library(Biostrings)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

# note this only makes a 10-mer context around 2bp deletions
muts[nchar(refnt) - nchar(altnt) == 2, deca := as.character(getSeq(BSgenome.Hsapiens.1000genomes.hs37d5, GRanges(seqnames=chr, ranges=IRanges(start=pos-3, end=pos+6))))]
muts[nchar(refnt) - nchar(altnt) == 2, deca.rc := ifelse(dn %in% will.revcomp, reverseComplement(DNAStringSet(deca)), DNAStringSet(deca))]

muts.tmp <- muts[nchar(refnt)-nchar(altnt)==2 & muttype=='indel' & age>50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & amp == 'PTA' & celltype=='neuron']

# Plot the first pink peak separately
#ggseqlogo(split(muts.tmp$deca.rc, paste(ifelse(muts.tmp$mutsig == '2:Del:R:0', 'Non-repeat', '1 or more repeat units'), muts.tmp$phenotype)), ncol=4) +
    #annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    #annotate('text', x=5.5, y=1.45, label='Deletion')

# Don't bother with the first pink peak, make 2x2 layout
ggseqlogo(split(muts.tmp[mutsig != '2:Del:R:0', deca.rc], muts.tmp[mutsig != '2:Del:R:0', phenotype]), ncol=2) +
    annotate('segment', x = 4.5, xend=6.5, y=1.3, yend=1.3, size=2) +
    annotate('text', x=5.5, y=1.45, label='Deletion') +
    labs(title='10-mer deletion context around 2 bp deletions', subtitle='First pink peak excluded, age>50, no MAPD/indel outliers') +
    xlab("Local deletion position")
dev.print(dev=pdf, file='10mer_deletion_context.NO_MAPD_INDEL_OUTLIERS_AGE_GT_50.pdf')


########################################################################################
#
# Code below copied from work/manual_enrichment/code.R
#
########################################################################################

enrich <- fread('manual_enrichment_analysis/enrichment.csv')

# quantile based covariates
q <- enrich[!is.na(QUANTILES) & BINSIZE=='1000']
q[, quantile := as.integer(quantile)]
q[datasource=='encode', datasource := ifelse(lineclass %in% c('H3K27me3', 'H3K9me3'), 'histone_inactive_mark', 'histone_active_mark')]
q[datasource=='gtex', datasource := 'gtex_expression']
q[datasource=='ucsc_track', datasource:= 'genetic conservation']
q[datasource== 'genetic conservation', datasource := paste('conservation', dataclass)]
q[, phenotype2 := paste(phenotype, delsig)]
# order the levels of delsig
q[, delsig := factor(delsig, levels=c('', 'negative', 'positive', 'positive_normal', 'positive_high'))]

q[, mapd.outlier := 'no_mapd_outliers']
# extra groups were only created when MAPD+indel outliers existed to remove
# loop through the added groups. for each group, there will be a corresponding group without
# the added fragment '_no_mapd_indel_outliers' that DOES contain outliers
for (gr in q[, unique(grep('_no_mapd_indel_outliers', group, value=T))]) {
    corr.grp <- sub(pattern='_no_mapd_indel_outliers', rep='', gr)
    q[group == corr.grp, mapd.outlier := 'with_mapd_outliers']
}

# filter down to just specific indel breakdowns
qq <- q[(datasource != 'gtex_expression' | dataclass=='Brain') & amp == 'PTA' & (delsig !='' & celltype == 'neuron') & (delsig %in% c('', 'negative', 'positive_normal','positive_high', 'positive'))]

# another filtration that does not break apart pink peak + and - cells
qqq <- q[(datasource != 'gtex_expression' | dataclass=='Brain') & amp == 'PTA' & delsig =='' & celltype == 'neuron']

#includes MDA
qqq2 <- q[(datasource != 'gtex_expression' | dataclass=='Brain') & delsig =='' & celltype == 'neuron']

# non-quantitative enrichments
b <- enrich[is.na(QUANTILES)]
b[, mapd.outlier := 'no_mapd_outliers']
# extra groups were only created when MAPD+indel outliers existed to remove
# loop through the added groups. for each group, there will be a corresponding group without
# the added fragment '_no_mapd_indel_outliers' that DOES contain outliers
for (gr in b[, unique(grep('_no_mapd_indel_outliers', group, value=T))]) {
    corr.grp <- sub(pattern='_no_mapd_indel_outliers', rep='', gr)
    b[group == corr.grp, mapd.outlier := 'with_mapd_outliers']
}
b[, mutsig.deltype := fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')]
########################################################################################
#
# END COPIED CODE
#
########################################################################################


#
# 11. Show the effect of the MAPD/indel outlier removal and the separate consideration of R0 pink peak
#
ggplot(qq[delsig=='positive_high' & !(muttype %in% c('A','indel_A')) & !grepl('_not_', muttype) & datasource==ds], aes(x=quantile, y=enr, col=phenotype, group=interaction(celltype, phenotype2, dataclass, lineclass))) +
    geom_line(alpha=0.1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    facet_grid(sub(pattern='indel_A_', rep='', muttype) ~ interaction(delsig, mapd.outlier), scales='free_y') +
    stat_summary(aes(group=interaction(celltype,phenotype2)), fun.y=mean, geom='line', linewidth=0.75) +
    base.theme +
    scale_y_log10() +
    labs(title="Quantitative enrichment - MAPD outlier and R0 pink peak effects", subtitle=ds) +
    ylab('Enrichment (obs/exp, log-scale)')



#
# 12. All enrichment analyses without previous bug (grep vs. grepl) in plotting, compare include MAPD vs. exclude MAPD outliers
#
for (ds in unique(qq$datasource)) {
    print(ds)
    p <- ggplot(qq[datasource==ds & !(muttype %in% c("A",'indel_A'))], aes(x=quantile, y=enr, col=phenotype, group=interaction(celltype, phenotype2, dataclass, lineclass))) +
        geom_line(alpha=0.1) +
        facet_grid(sub(pattern='indel_A_', rep='', muttype) ~ forcats::fct_cross(mapd.outlier, delsig, sep='\n'), scales='free_y') +
        stat_summary(aes(group=interaction(celltype,phenotype2)), fun.y=mean, geom='line', linewidth=0.75) +
        base.theme +
        geom_hline(yintercept=1, linewidth=0.25) +
        scale_y_log10() +
        labs(title="Quantitative enrichment", subtitle=ds) +
        ylab('Enrichment (obs/exp, log-scale)')
    ggsave(filename=paste0('manual_enrichment_analysis/manual_plots/', ds, '.pdf'), plot=p, height=16, width=11)
}



#
# 13. Same as above, but do not show MAPD+indel outliers
#
for (ds in unique(qq$datasource)) {
    print(ds)
    p <- ggplot(qq[mapd.outlier=='no_mapd_outliers' & datasource==ds & !(muttype %in% c("A",'indel_A'))], aes(x=quantile, y=enr, col=phenotype, group=interaction(celltype, phenotype2, dataclass, lineclass))) +
        geom_line(alpha=0.1) +
        facet_grid(sub(pattern='indel_A_', rep='', muttype) ~ forcats::fct_cross(mapd.outlier, delsig, sep='\n'), scales='free_y') +
        stat_summary(aes(group=interaction(celltype,phenotype2)), fun.y=mean, geom='line', linewidth=0.75) +
        base.theme +
        geom_hline(yintercept=1, linewidth=0.25) +
        scale_y_log10() +
        labs(title="Quantitative enrichment", subtitle=ds) +
        ylab('Enrichment (obs/exp, log-scale)')
    ggsave(filename=paste0('manual_enrichment_analysis/manual_plots/', ds, '.NO_MAPD_INDEL_OUTLIERS.pdf'), plot=p, height=16, width=11)
}


#
# 13a. MAIN FIGURE panels: 2x2 table of enrichment vs. GTEx expression
#
ggplot(qq[mapd.outlier=='no_mapd_outliers' & datasource=='gtex_expression' & delsig %in% c('negative', 'positive') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), group=interaction(celltype, phenotype2, dataclass, lineclass))) +
    geom_line(alpha=0.1) +
    facet_grid(fct_shift(fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')) ~ fct_recode(delsig, `ID83A positive`='positive', `ID83A negative`='negative')) +
    stat_summary(aes(group=interaction(celltype,phenotype2)), fun.y=mean, geom='line', linewidth=0.75) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes vs. GTEx expression", subtitle='Thin lines: individual brain tissues, thick lines: cross-brain tissue average') +
    xlab("Gene expression quantile") + ylab("Indel enrichment (obs / exp, log-scale)")
dev.print(dev=pdf, file='enrichment_gtex_expression_2x2.pdf')

#
# 13b. MAIN FIGURE panels: 2x2 table of enrichment vs. GTEx expression, but no per-tissue thin lines
#      to improve readability.
#
ggplot(qq[mapd.outlier=='no_mapd_outliers' & datasource=='gtex_expression' & delsig %in% c('negative', 'positive') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), group=interaction(celltype, phenotype2, dataclass, lineclass))) +
    facet_grid(fct_shift(fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')) ~ fct_recode(delsig, `ID83A positive`='positive', `ID83A negative`='negative')) +
    stat_summary(aes(group=interaction(celltype,phenotype2)), fun.y=mean, geom='line', linewidth=0.75) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes vs. GTEx expression", subtitle='Only cross-brain-tissue averages to emphasize positive trends') +
    xlab("Gene expression quantile") + ylab("Indel enrichment (obs / exp, log-scale)")
dev.print(dev=pdf, file='enrichment_gtex_expression_2x2.NO_TISSUE_SPECIFIC_LINES.pdf')


#
# 13c. MAIN FIGURE panels: 2x2 table of enrichment vs. GTEx expression, but using only the PFC (BA9)
#      GTEx expression data. All of our neurons are either BA9 or BA6; previous papers used only BA9,
#      so using this one tissue was more justifiable.
#
ggplot(qq[mapd.outlier=='no_mapd_outliers' & datasource=='gtex_expression' & lineclass=='Brain_-_Frontal_Cortex__BA9_' & delsig %in% c('negative', 'positive') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), group=interaction(celltype, phenotype2, dataclass, lineclass))) +
    geom_line(linewidth=0.75) +
    facet_grid(fct_shift(fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')) ~ fct_recode(delsig, `ID83A positive`='positive', `ID83A negative`='negative')) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes vs. GTEx expression", subtitle='Only PFC (BA9) GTEx expression data') +
    xlab("Gene expression quantile") + ylab("Indel enrichment (obs / exp, log-scale)")
dev.print(dev=pdf, file='enrichment_gtex_expression_2x2.JUST_BA9_LINE.pdf')


#
# 13d. MAIN FIGURE panels: enrichment vs. GTEx expression, but using only the PFC (BA9)
#      GTEx expression data. This time, drop the panel for disease-associated deletions in
#      ID83A-negative cells, where due to the small counts the lines are just noise.
#
ggplot(qq[mapd.outlier=='no_mapd_outliers' & datasource=='gtex_expression' & lineclass=='Brain_-_Frontal_Cortex__BA9_' & delsig %in% c('negative', 'positive') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control') & !(delsig=='negative' & muttype=='indel_A_2del_3delMH')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'))) +
    geom_line(linewidth=0.75) +
    facet_wrap(~fct_cross(fct_shift(fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')), fct_recode(delsig, `ID83A positive`='positive', `ID83A negative`='negative'), sep='\n'), ncol=3) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes vs. GTEx expression", subtitle='Only PFC (BA9) GTEx expression data') +
    xlab("Gene expression quantile") + ylab("Indel enrichment (obs / exp, log-scale)")
dev.print(dev=pdf, file='enrichment_gtex_expression_2x2.JUST_BA9_LINE.3_PANELS_ONLY.pdf')

#
# 13e.2. MAIN FIGURE panels: enrichment vs. GTEx expression, but using only the PFC (BA9)
#      GTEx expression data. No longer separate by pink peak +/- status because this reduces
#      N too much.
#
#      Now just comparing aging-related indels vs. disease-related deletions for each phenotype.
#      One panel per phenotype. Now that there are so few lines, it is feasible to add signif
#      asterisks, which reveal that even when disease-related deletions appear to positively
#      correlate with expression, the enrichments are usually not significant (or less signif
#      than aging indels for the same quantile).
#
#   e.2 UPDATE: now show the scRNAseq for excitatory neurons rather than GTEx in the main figure.
#      supplement now contains GTEx.
#
qqq[, mutsig.deltype := fct_recode(muttype, `Aging-related indels`='indel_A_not_2del_3delMH', `Disease-related deletions`='indel_A_2del_3delMH')]
qqq[, phenotype := factor(phenotype, levels=c('control', 'als', 'ftd', 'alz'))] # order matters here
# this version (with fct_shift()) gives a solid line for aging, dotted for diseased
#ggplot(qqq[mapd.outlier=='no_mapd_outliers' & datasource=='gtex_expression' & lineclass=='Brain_-_Frontal_Cortex__BA9_' & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control') & !(delsig=='negative' & muttype=='indel_A_2del_3delMH')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), linetype=fct_shift(fct_drop(mutsig.deltype)))) +
ggplot(qqq[mapd.outlier=='no_mapd_outliers' & datasource=='scrnaseq' & lineclass=='Excitatory-Neurons' & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & !(delsig=='positive' & phenotype=='control') & !(delsig=='negative' & muttype=='indel_A_2del_3delMH')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), linetype=fct_drop(mutsig.deltype))) +
    geom_line(linewidth=0.5) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    geom_point(aes(shape=ifelse(pval < 0.01, 'P < 0.01', 'n.s.')), fill='white', stroke=2/4, size=1.5) +
    facet_wrap(~fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), ncol=4) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes vs. snRNA-seq expression", subtitle='Only excitatory neurons from Ganz et al., * < 0.01, ** < 0.001, *** < 0.0001', shape='') +
    xlab("Gene expression decile") + ylab("Enrichment\n(obs / exp, log-scale)") +
    scale_shape_manual(values=c(`P < 0.01`=16, `n.s.`=21), breaks=c('P < 0.01', 'n.s.')) +
    scale_x_continuous(breaks=1:10)

    #geom_text_repel(aes(label=sapply(pmax(0, floor(-log10(pval))-1), function(i) paste(rep('*', i), collapse=''))), segment.size=0.25, min.segment.length=0, family='mono') #, box.padding=0.3)
    #geom_text(aes(label=sapply(pmax(0, floor(-log10(pval))-1), function(i) paste(rep('*', i), collapse=''))), nudge_y=1/50)
dev.print(dev=pdf, file='enrichment_scrnaseq.JUST_EXC_NEU.NO_SEPARATION_OF_POS_AND_NEG.SIGNIF_NOT_BY_STARS.SOLID_LINE_FOR_DISEASE.pdf')



#
# 13f.2. SUPPLEMENTARY FIGURE panels: enrichment using only the PFC (BA9), not separated by
#      pink peak+ or - status.
#
#      THIN LINES and THICK LINES version - thus no p-values
#
#   f.2 UPDATE: now show GTEx for BA9 in the supplement, scRNAseq in the main figure.
#
ggplot(qqq[mapd.outlier=='no_mapd_outliers' & !(datasource %in% c('repliseq', 'scrnaseq', 'conservation phastCons', 'conservation phyloP100way')) & (datasource != 'gtex_expression' | lineclass=='Brain_-_Frontal_Cortex__BA9_') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), linetype=fct_drop(mutsig.deltype), group=interaction(fct_drop(mutsig.deltype), lineclass))) +
    geom_line(alpha=0.1, linewidth=0.5) +
    stat_summary(aes(group=interaction(fct_drop(mutsig.deltype))), fun.y=mean, geom='line', linewidth=0.75) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    facet_grid(datasource~fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz')) +
    scale_y_log10() +
    labs(title="Quantitative enrichment: indel subtypes", subtitle='Thin lines: separate marks or cell types; thick lines: average of thin lines') +
    xlab("Covariate quantile (higher = more transcriptional activity)") + ylab("Indel enrichment (obs / exp, log-scale)")
dev.print(dev=pdf, file='enrichment_supplement_quantitative.pdf')



#
# 13g.2. SUPPLEMENTARY FIGURE panels: not separated by pink peak+ or - status.
#
#      Only show exc. neu. lines for scATAC/RNA-seq; separate rows for histone marks
#
#   g.2 UPDATE: swap GTEx and scrnaseq/exc.neu. in main and supp figs.
#
ggplot(qqq[mapd.outlier=='no_mapd_outliers' & !(datasource %in% c('repliseq', 'scrnaseq', 'conservation phastCons', 'conservation phyloP100way')) & (datasource != 'gtex_expression' | lineclass=='Brain_-_Frontal_Cortex__BA9_') & (datasource != 'scatacseq' | lineclass == 'excitatory_neuron') & (datasource != 'histone_inactive_mark' | lineclass == 'H3K9me3') & (datasource != 'histone_active_mark' | lineclass != 'H3K36me3') & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH')], aes(x=quantile, y=enr, col=fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz'), linetype=fct_drop(mutsig.deltype), group=interaction(fct_drop(mutsig.deltype), lineclass))) +
    geom_line(linewidth=0.3) +
    geom_point(aes(shape=ifelse(pval < 0.01, 'P < 0.01', 'n.s.')), fill='white', stroke=2/4, size=1/2) +
    ggtheme + theme(aspect.ratio=1) +
    geom_hline(yintercept=1, linewidth=0.25) +
    facet_grid(fct_cross(datasource,lineclass,sep='\n')~fct_recode(phenotype, Control='control', ALS='als', FTD='ftd', Alzheimers='alz')) +
    scale_y_log10() +
    xlab("Covariate decile") + ylab("Indel enrichment (obs / exp, log-scale)") +
    theme(strip.text.y.right=element_text(angle=0), axis.text.x=element_text(size=8)) +
    scale_x_continuous(breaks=1:10) +
    scale_shape_manual(values=c(`P < 0.01`=16, `n.s.`=21), breaks=c('P < 0.01', 'n.s.'))

    #geom_text_repel(aes(label=sapply(pmax(0, floor(-log10(pval))-1), function(i) paste(rep('*', i), collapse=''))), segment.size=0.15, min.segment.length=0, size=2.5, max.overlaps=13, family='mono') + 
    #labs(title="Quantitative enrichment: indel subtypes") +
    #stat_summary(aes(group=interaction(fct_drop(mutsig.deltype))), fun.y=mean, geom='line', linewidth=0.75) +
    #geom_text(aes(label=sapply(pmax(0, floor(-log10(pval))-1), function(i) paste(rep('*', i), collapse=''))), nudge_y=1/50) +
dev.print(dev=pdf, file='enrichment_supplement_quantitative.SELECTED_TRACKS_WITH_PVALUES.NO_ASTERISKS_FOR_PVALUES.pdf')


#
# 13h. SUPPLEMENTARY FIGURE panels: NON-QUANTITATIVE analyses
#
# N.B.: I removed non-neuronal enh/prom data for brevity. ChromHMM states were removed for:
#   (1) being too small (e.g., 10_TssBiv, ...) or (2) having no significant p-value across
#   all phenotypes (e.g., 2_TssAFlnk..)
bb <- b[datasource %in% c('gtex', 'nott_enhprom', 'roadmap') & !(quantile %in% c('10_TssBiv', '11_BivFlnk', '12_EnhBiv', '3_TxFlnk', '8_ZNF/Rpts', '13_ReprPC', '2_TssAFlnk', '4_Tx', '6_EnhG', 'outside')) & (datasource != 'nott_enhprom' | quantile %in% c('neuron_promoter', 'neuron_enhancer')) & mapd.outlier=='no_mapd_outliers' & muttype %in% c('indel_A_2del_3delMH', 'indel_A_not_2del_3delMH') & delsig=='' & amp=='PTA' & celltype=='neuron']

bb[, quantile := fct_recode(quantile, Transcribed='tx', Untranscribed='utx')]
bb[, quantile := fct_relevel(quantile, c('1_TssA', '5_TxWk', '7_Enh', '9_Het', '14_ReprPCWk', '15_Quies'))]
bb[, datasource := ifelse(datasource != 'roadmap', datasource, ifelse(quantile %in% c('1_TssA', '5_TxWk', '7_Enh'), 'ChromHMM active', 'ChromHMM inactive'))]
bb[datasource == 'gtex', datasource := 'Gene']
bb[datasource == 'nott_enhprom', datasource := 'Regulatory element']

ggplot(bb, aes(x=quantile, y=enr, fill=mutsig.deltype)) +
    geom_col(position='dodge') +
    geom_hline(yintercept=1, linewidth=0.25) +
    facet_grid(phenotype~datasource, scale='free_x', space='free_x', labeller=label_wrap_gen(width=12)) +
    scale_y_log10() +
    base.theme +
    geom_text(aes(label=sapply(pmax(0, floor(-log10(pval))-1), function(i) paste(rep('*', i), collapse=''))), position=position_dodge(width=0.9)) +
    theme(axis.text.x=element_text(angle=90, vjust=1/2, hjust=1)) +
    xlab(element_blank()) + ylab('Indel enrichment (obs/exp, log-scale)')
dev.print(dev=pdf, file='enrichment_supplement_bed_regions.SELECTED_TRACKS_WITH_PVALUES.pdf')


#
# 14. Indel size distribution plot
#
muts <- meta[fread('tables/all___FILTERED_mut___any.csv'),, on=.(sample)]

ggplot(muts[muttype == 'indel' & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA' & abs(nchar(refnt)-nchar(altnt))<20], aes(x=nchar(altnt)-nchar(refnt), col=phenotype, fill=phenotype)) +
    geom_bar(width=1/20) +
    facet_wrap(~phenotype, ncol=2, scale='free_y') +
    geom_vline(xintercept=0, linewidth=0.2, linetype=8) +
    ylab('Number of sIndels') + xlab('Indel size') +
    labs(title='Indel sizes', subtitle='PTA neurons only, |size|<20 bp') +
    ggtheme + theme(aspect.ratio=2/3)
dev.print(dev=pdf, file='indel_size_distribution.NOT_RIDGE_PLOT.pdf')



#
# 15. De novo SBS signatures
#
sbs.sigs <- fread('mutsigs/sigprofilerextractor/all/SBS96/Most_Stab_Sigs/signatures.csv')
# convert COSMIC format to SCAN2 matrix format
sbs.sigs.mat <- as.matrix(sbs.sigs[,-1],
    rownames=sbs.sigs[, paste0(substr(MutationType,1,1), substr(MutationType,3,3), substr(MutationType,7,7), ':', substr(MutationType,3,5))])[names(table(sbs96(c()))),]
plot.sbs96(sbs.sigs.mat, show.detailed.types=T, uniform.y.axis=T)
dev.print(dev=pdf, file='de_novo_sbs96_signatures.pdf')



#
# 16. De novo ID signatures
#
id.sigs <- fread('mutsigs/sigprofilerextractor/all/ID83_corrected/Most_Stab_Sigs/signatures.csv')
# convert COSMIC format to SCAN2 matrix format
id.sigs.mat <- as.matrix(id.sigs, rownames=1)[names(table(id83(c()))),]
# cbind(., 0) is just a trick to get the SBS and ID plots to look similar. R's layout
# changes text/line sizes with there are >=3 panels in a plot.
plot.id83(cbind(id.sigs.mat,0), show.detailed.types=F, uniform.y.axis=F, las=2)
dev.print(dev=pdf, file='de_novo_id83_signatures.pdf')



#
# 17. Spectra by cell type/phenotype and pink peak class
#
muts2 <- muts[amp == 'PTA' & muttype == 'indel' & sample %in% no.mapd.indel.outliers]
muts2.split <- split(muts2, paste(ifelse(muts2$delsig=='negative', 'negative', 'positive'), muts2$phenotype, muts2$celltype))
muts2.split <- muts2.split[c(1:3, 5, 4, 6:9)]  # hack to reorder the control oligo panel
plot.id83(sapply(muts2.split, function(df) table(id83(df$mutsig))), max.nrow=5, uniform.y.axis=F)
dev.print(dev=pdf, file='id83_spectra_neg_vs_pos.pdf')



########################################################################
#
# Outliers, statistics, tables
#
########################################################################

# These false positive rates were derived in Luquette et all, Nat Genet 2022
# from synthetic data with known spike-in somatic mutations at various rates.
# IMPORTANT: these FPRs were estimated without the post-processing, cohort-wide
# step of recurrence filtering, which likely reduces false positives.  Thus,
# all estimated false positive numbers are likely overestimates.
#
# FPRs in false calls per million analyzable basepairs.
#
# These FPRs were also used in Ganz, Luquette et al, Cell 2024 to estimate
# per-cell and total false discovery rates.
fpr.per.mb.snv=0.01313788
fpr.per.mb.indel=0.0007298824

summary.objects <- list.files(pattern='.*.rda', path='scan2/summary_objects', full.name=T)

if (reload) {
    per.cell.stats <- rbindlist(lapply(summary.objects, function(fn) {
        print(fn)
        load(fn)
        mapd <- mapd(smry, type='canonical')

        tdata <- training(smry)
        # Not the entry in the mutburden data frame. That entry estimates the sensitivity of
        # the genome with moderate sequencing depth (middle 50% of the distribution).
        snv.sens <- tdata[resampled.training.site==T & muttype=='snv', mean(training.pass)]
        indel.sens <- tdata[resampled.training.site==T & muttype=='indel', mean(training.pass)]

        dptab <- decompress.dt(smry@depth.profile$dptab)
        mean.depth <- sum(rowSums(dptab) * as.numeric(rownames(dptab))) / sum(dptab)
        # WARNING! This is the SNV minimum - indels are slightly higher, so the fraction of
        # genome analyzable for indels is smaller.
        sfp <- smry@static.filters$params$snv
        snv.analyzable.basepairs <- sum(dptab[(sfp$min.sc.dp+1):nrow(dptab), (sfp$min.bulk.dp+1):ncol(dptab)])
        sfp <- smry@static.filters$params$indel
        indel.analyzable.basepairs <- sum(dptab[(sfp$min.sc.dp+1):nrow(dptab), (sfp$min.bulk.dp+1):ncol(dptab)])

        pdata <- passing(smry, passtype='vafbased')
        snv.calls <- pdata[muttype=='snv', sum(pass)]
        snv.filtered.calls <- muts[muttype=='snv' & sample==name(smry), sum(pass)]
        snv.fps <- min(snv.filtered.calls, snv.analyzable.basepairs/1e6 * fpr.per.mb.snv)
        snv.fdr <- ifelse(snv.filtered.calls == 0, 0, snv.fps/snv.filtered.calls)
        indel.calls <- pdata[muttype=='indel', sum(pass)]
        indel.filtered.calls <- muts[muttype=='indel' & sample==name(smry), sum(pass)]
        indel.fps <- min(indel.filtered.calls, indel.analyzable.basepairs/1e6 * fpr.per.mb.indel)
        indel.fdr <- ifelse(indel.filtered.calls == 0, 0, indel.fps/indel.filtered.calls)

        rbind(
            # not muttype dependent
            melt(data.table(sample=name(smry), muttype=NA, mapd=mapd, sequencing.depth=mean.depth), id.vars=1:2),
            # fraction.genome.analyzable - divide by 2 because get.gbp.by.genome returns sum of gbp
            #   over maternal and paternal copies.
            melt(data.table(sample=name(smry), muttype='snv', sens=snv.sens, analyzable.basepairs=snv.analyzable.basepairs, fraction.genome.analyzable=snv.analyzable.basepairs/(1e9*get.gbp.by.genome(smry)/2), calls=snv.calls, filtered.calls=snv.filtered.calls, fps=snv.fps, fdr=snv.fdr), id.vars=1:2),
            melt(data.table(sample=name(smry), muttype='indel', sens=indel.sens, analyzable.basepairs=indel.analyzable.basepairs, fraction.genome.analyzable=indel.analyzable.basepairs/(1e9*get.gbp.by.genome(smry)/2), calls=indel.calls, filtered.calls=indel.filtered.calls, fps=indel.fps, fdr=indel.fdr), id.vars=1:2)
        )
    }))
    per.cell.stats <- meta[,.(sample,donor,sex,age,phenotype,amp,tissue.origin,tdp43,celltype,`ID-A status`=fct_collapse(delsig, positive=c('normal', 'high')))][per.cell.stats,,on=.(sample)]
    per.cell.stats[, group := fct_cross(celltype, phenotype, amp, sep=' ')]
    per.cell.stats[, muttype := fct_shift(fct_recode(muttype, SNV='snv', Indel='indel'))]
}


#
# 18. Aggregate stats
#
p1 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & is.na(muttype)], aes(x=group, y=value, fill=phenotype)) + facet_wrap(~fct_recode(variable, MAPD='mapd', `Seq. depth`='sequencing.depth'), scale='free_y', switch='y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/5, width=0.2) + xlab(element_blank()) + ylab(element_blank()) + theme(strip.placement='outside')
p2 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & variable == 'filtered.calls'], aes(x=group, weight=value, fill=phenotype)) + facet_wrap(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_bar() + xlab(element_blank()) + ylab('Mutation calls')
(p1|p2)+plot_layout(guides = "collect")&theme(aspect.ratio=1)
dev.print(dev=pdf, file='per_cell_stats.global.WITH_MDA.pdf')


#
# 18b. Aggregate stats, MDA columns excluded
#
p1 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & amp == 'PTA' & is.na(muttype)], aes(x=group, y=value, fill=phenotype)) + facet_wrap(~fct_recode(variable, MAPD='mapd', `Seq. depth`='sequencing.depth'), scale='free_y', switch='y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/5, width=0.2) + xlab(element_blank()) + ylab(element_blank()) + theme(strip.placement='outside')
p2 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable == 'filtered.calls'], aes(x=group, weight=value, fill=phenotype)) + facet_wrap(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_bar() + xlab(element_blank()) + ylab('Mutation calls')
(p1|p2)+plot_layout(guides = "collect") & theme(aspect.ratio=1)
dev.print(dev=pdf, file='per_cell_stats.global.pdf')

#
# 18c. Aggregate stats, MDA columns excluded, age matched > 50
#
p1 <- ggplot(per.cell.stats[age > 50 & sample %in% no.mapd.indel.outliers & amp == 'PTA' & is.na(muttype)], aes(x=group, y=value, fill=phenotype)) + facet_wrap(~fct_recode(variable, MAPD='mapd', `Seq. depth`='sequencing.depth'), scale='free_y', switch='y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/5, width=0.2) + xlab(element_blank()) + ylab(element_blank()) + theme(strip.placement='outside')
p2 <- ggplot(per.cell.stats[age > 50 & sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable == 'filtered.calls'], aes(x=group, weight=value, fill=phenotype)) + facet_wrap(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_bar() + xlab(element_blank()) + ylab('Mutation calls')
(p1|p2)+plot_layout(guides = "collect") & theme(aspect.ratio=1)
dev.print(dev=pdf, file='per_cell_stats.global.AGE_MATCHED_DO_NOT_USE_BECAUSE_NEW_BA6_NOT_INCLUDED.pdf')


#
# 19a. Per-cell stats, MDA columns excluded, all ages
#
p1 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable %in% c('fraction.genome.analyzable', 'sens')], aes(x=group, y=value, fill=phenotype)) + facet_grid(fct_recode(variable, Sensitivity='sens', `Fraction genome\n>=min seq. depth`='fraction.genome.analyzable') ~ muttype, scale='free_y', switch='y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/5, width=0.2) + xlab(element_blank()) + theme(strip.placement='outside') + ylab(element_blank())

fp.tab <- dcast(per.cell.stats[sample %in% no.mapd.indel.outliers & variable %in% c('fps', 'filtered.calls')], amp + phenotype + celltype + group + muttype ~ variable, value.var='value', fun.aggregate=sum)
fp.tab[, catalog.fdr := fps/filtered.calls]
fp.tab <- melt(fp.tab, id.vars=1:5, measure.vars='catalog.fdr')

p2 <- ggplot(per.cell.stats[sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable=='fdr'], aes(x=group, y=value, fill=phenotype)) + facet_grid(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/5, width=0.2) + ylab("False discovery rate")

# does not need no.mapd.indel.outliers because fp.tab was already built with it
p3 <- ggplot(fp.tab[amp == 'PTA' & variable=='catalog.fdr'], aes(x=group, weight=value, fill=phenotype)) + facet_grid(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_bar() + xlab(element_blank()) + geom_text(stat='count', aes(label=after_stat(round(count,2))), nudge_y=0.006, size=3) + ylab('Catalog-wide FDR')

(free(p1+guides(fill='none')) + (p2/p3 + plot_layout(axes='collect')))
dev.print(dev=pdf, file='per_cell_stats.pdf')


#
# 19b. Per-cell stats, MDA columns excluded, age matched using age>50
#
p1 <- ggplot(per.cell.stats[age > 50 & sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable %in% c('fraction.genome.analyzable', 'sens')], aes(x=group, y=value, fill=phenotype)) + facet_grid(fct_recode(variable, Sensitivity='sens', `Fraction genome\n>=min seq. depth`='fraction.genome.analyzable') ~ muttype, scale='free_y', switch='y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/2, width=0.2) + xlab(element_blank()) + theme(strip.placement='outside') + ylab(element_blank())

fp.tab <- dcast(per.cell.stats[age > 50 & sample %in% no.mapd.indel.outliers & variable %in% c('fps', 'filtered.calls')], amp + phenotype + celltype + group + muttype ~ variable, value.var='value', fun.aggregate=sum)
fp.tab[, catalog.fdr := fps/filtered.calls]
fp.tab <- melt(fp.tab, id.vars=1:5, measure.vars='catalog.fdr')

p2 <- ggplot(per.cell.stats[age > 50 & sample %in% no.mapd.indel.outliers & amp == 'PTA' & variable=='fdr'], aes(x=group, y=value, fill=phenotype)) + facet_grid(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=1/2, width=0.2) + ylab("False discovery rate")#xlab(element_blank()) + theme(strip.placement='outside') + ylab(element_blank())

p3 <- ggplot(fp.tab[age > 50 & amp == 'PTA' & variable=='catalog.fdr'], aes(x=group, weight=value, fill=phenotype)) + facet_grid(~muttype, scale='free_y') + ggtheme + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2)) + geom_bar() + xlab(element_blank()) + geom_text(stat='count', aes(label=after_stat(round(count,2))), nudge_y=0.006) + ylab('Catalog-wide FDR')

(free(p1+guides(fill='none')) + (p2/p3 + plot_layout(axes='collect')))
dev.print(dev=pdf, file='per_cell_stats.AGE_MATCHED.pdf')


per.cell.stats[variable=='mapd',.(sample,mapd=value)][resids[burden.source=='denovo' & burden.type=='SBS96' & burden.subtype=='SBS96B', .(sample, burden)],,on=.(sample)]




#
# 21a. Quality check - are SBS-B and/or ID-A artifactual? One possible indicator would be
#   strong correlation with a quality metric, e.g., MAPD.
#

qc.tab <- per.cell.stats[variable=='mapd',.(sample,mapd=value)][resids[burden.source=='denovo' & ((burden.type=='SBS96' & burden.subtype=='SBS96B') | (burden.type=='ID83' & burden.subtype=='ID83A')), .(donor, sample, age, amp, phenotype, delsig, celltype, burden.subtype, burden, expected.burden)][sample %in% no.mapd.indel.outliers & age > 50],,on=.(sample)]

ggplot(qc.tab, aes(x=mapd, y=burden/expected.burden, col=phenotype)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE, linewidth=2/3, col='black') +
    ggpmisc::stat_fit_glance(method='lm', aes(label=sprintf('P = %0.2g, R^2 = %0.2g', after_stat(p.value), after_stat(r.squared))), size=2.5) +
    facet_grid(burden.subtype~phenotype+amp+celltype, scale='free', switch='y') +
    theme(aspect.ratio=1, strip.placement='outside') +
    ylab('Obs/exp') 
dev.print(dev=pdf, file='mapd_vs_sbsB_and_id83A.pdf')

#
# 21b. Split by deletion signature, colors according to donor
#
ggplot(qc.tab, aes(x=mapd, y=burden/expected.burden, col=phenotype)) +
    geom_point(aes(col=donor)) +
    geom_smooth(method='lm', se=FALSE, linewidth=2/3, col='black') +
    ggpmisc::stat_fit_glance(method='lm', aes(label=sprintf('P = %0.2g, R^2 = %0.2g', after_stat(p.value), after_stat(r.squared))), size=2.5) +
    facet_grid(burden.subtype+delsig~phenotype+amp+celltype, scale='free', switch='y') +
    theme(aspect.ratio=1, strip.placement='outside') +
    ylab('Obs/exp') 
dev.print(dev=pdf, file='mapd_vs_sbsB_and_id83A.BY_DELSIG_AND_DONOR.pdf')



disease.related <- c(
     "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",  # pink peaks
     "2:Del:M:1", "3:Del:M:1", "3:Del:M:2"  # microhomology purple peaks
)



#
# 22. COSMIC ID4 - just plotting ID4 from COSMIC
#
plot.id83(as.matrix(fread('v3.3_ID4_PROFILE.txt'), rownames=1))
dev.print(dev=pdf, file='cosmic_id4.pdf')


#
# 23. Stats for catalog size, mean sens and FDR
#

# all mut count - includes OLs
muts[rescue == FALSE & sample %in% no.mapd.indel.outliers & amp=='PTA', table(muttype)]



#
# 24. COSMIC SBS30 vs. ID4 comparison
#
#   To sidestep any issues about SBS-B only surfacing from de novo extraction
#   of PTA+MDA cells, just do direct COSMIC fitting with SigProfilerAssignment.
#
spa <- rbind(
    melt(fread('manual_sigprofiler_assignment/all/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'), id.vars=1),
    melt(fread('manual_sigprofiler_assignment/all_ID83_corrected//Assignment_Solution/Activities/Assignment_Solution_Activities.txt'), id.vars=1)
)
colnames(spa)[1] <- 'sample'
spa <- meta[spa,,on=.(sample)]
spa[sample %in% no.mapd.indel.outliers, sig.lump := fct_lump_min(variable, min=200, w=value)]
# Shows ID4+1 vs SBS30
ggplot(dcast(spa[age > 50 & celltype=='neuron' & amp =='PTA' & sample %in% no.mapd.indel.outliers], sample+amp+celltype+phenotype+age ~ variable, value.var='value'), aes(x=ID4+1, y=SBS30, col=phenotype)) + geom_point() + geom_smooth(method='lm', se=F) + scale_x_log10()


# 25. Fraction of 2bp dels among all indels
muts.new.ol.phenotype <- copy(muts)[, phenotype := fct_relevel(factor(ifelse(celltype=='oligo', 'Control oligo', ifelse(celltype=='neuron' & phenotype=='Control', 'Control neuron', as.character(phenotype)))), c('Control neuron', 'Control oligo', 'ALS', 'FTD', 'Alzheimers'))]
fig3b <- ggplot(muts.new.ol.phenotype[age > 50 & muttype == 'indel' & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & amp == 'PTA', .(frac.2bp=mean(nchar(refnt)-nchar(altnt)==2)), by=.(sample, phenotype)], aes(x=fct_drop(phenotype), y=frac.2bp, fill=fct_drop(phenotype))) +
    geom_boxplot(outliers=FALSE) +
    geom_jitter(size=1/4, width=0.2) +
    ylab('#2 bp deletions / all sIndels') + xlab(element_blank()) +
    labs(title='Indel sizes', subtitle='PTA neurons only, age>50, no MAPD/indel outliers') +
    ggtheme + theme(aspect.ratio=5/4) +
    expand_limits(y=0:1) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    geom_signif(comparisons=list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')), margin_top=-0.35, textsize=2.5, step_increase=1/10, size=1/4)
dev.print(dev=pdf, file='indel_2_bp_deletion_fraction_and_stats.pdf')


#
# 26. Fraction of pink peak positive cells with odds ratio test
# 
class.tab2 <- meta[amp=='PTA' & sample %in% no.mapd.indel.outliers,.(`ID-A high`=sum(delsig=='high'|delsig=='normal'), `ID-A normal`=sum(delsig=='negative'), total=nrow(.SD), rate=mean(delsig!='negative')),by=.(phenotype=factor(ifelse(phenotype=='Control',ifelse(celltype=='oligo', 'Control oligo', 'Control neuron'), as.character(phenotype)), levels=c('Control neuron', 'Control oligo', 'ALS', 'FTD', 'Alzheimers')))]

fishertests2 <- rbindlist(lapply(list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')), function(comp) {
    m <- as.matrix(class.tab2[phenotype %in% comp, .(`ID-A normal`, `ID-A high`)])  # normal first
    ft <- fisher.test(m)
    ci <- ft$conf.int
    data.frame(pheno1=comp[1], pheno2=comp[2],
        pheno1.n=sum(m[1,]), pheno1.positive=m[1,'ID-A high'],
        pheno2.n=sum(m[2,]), pheno2.positive=m[2,'ID-A high'],
        p.value=ft$p.value, odds.ratio=ft$estimate, ci.lower=ci[1], ci.upper=ci[2], ci.level=attr(ci, 'conf.level'))        
}))
fishertests2[, pheno1.frac := pheno1.positive/pheno1.n]
fishertests2[, pheno2.frac := pheno2.positive/pheno2.n]

# R's fisher.test() will not allow adding 0.5, so compute it manually
# So the OR here will not exactly match the CI computed by fisher.test
fishertests2[, odds.ratio.haldane.anscombe := ((pheno2.positive+0.5)/(pheno2.n-pheno2.positive+0.5))/((pheno1.positive+0.5)/(pheno1.n-pheno1.positive+0.5))]

fishertests2[, pheno2 := factor(pheno2, levels=c('Control oligo', 'ALS', 'FTD', 'Alzheimers'))]

ggplot(class.tab2, aes(x=phenotype, weight=rate, fill=phenotype, label=paste(`ID-A high`, total, sep='/'))) +
    geom_bar() +
    geom_text(aes(y=rate), nudge_y=0.035, size=3) +
    ggtheme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), aspect.ratio=5/4) +
    labs(title='Rates of ID-A high neurons', subtitle='PTA only, all ages, no MAPD/indel outliers') +
    xlab(element_blank()) + ylab('Fraction ID-A high cells') +
    expand_limits(y=c(0,1)) +
    geom_signif(aes(y=rate), comparisons=list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')),
        step_increase=1/12, margin_top=-0.21, textsize=2.75, size=0.15, annotation=sprintf('%0.2g', fishertests2$p.value), position=position_nudge(y=-0.05))
dev.print(dev=pdf, file='rate_of_id4_positive_cells_per_phenotype.pdf')


#
# 27. Odds ratio direct plot (i.e., not odds ratio P-values on %2bp del plot)
#
ggplot(fishertests2, aes(y=fct_rev(pheno2), x=odds.ratio.haldane.anscombe, label=round(odds.ratio.haldane.anscombe,1), col=pheno2)) +
    geom_point(size=2.5) +
    geom_errorbarh(aes(xmin=ci.lower, xmax=ci.upper), height=0.2) +
    geom_text(nudge_y=0.35, size=3, col='black') +
    ggtheme +
    xlab("Odds ratio vs. control neurons") +
    scale_x_log10() +
    theme(aspect.ratio=10/16) +
    geom_vline(xintercept=1, linetype=2) +
    ylab(element_blank()) +
    labs(title='Odds ratios and 95% conf ints', subtitle='Haldane-Anscombe correction (+0.5) to deal\nwith 0 oligo counts, no MAPD/indel outliers')
dev.print(dev=pdf, file='rate_of_id4_positive_cells_per_phenotype.ODDS_RATIO.pdf')




#############################################################################
#
# Figure 3 patchwork attempt
#
#############################################################################

# A. Indel size distribution plot
fig3a <- ggplot(muts[muttype == 'indel' & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA' & abs(nchar(refnt)-nchar(altnt))<20], aes(x=nchar(altnt)-nchar(refnt), col=phenotype, fill=phenotype)) +
    geom_bar(width=1/20) +
    facet_wrap(~phenotype, ncol=2, scale='free_y') +
    geom_vline(xintercept=0, linewidth=0.2, linetype=8) +
    ylab('Number of sIndels') + xlab('Indel size') +
    labs(title='Indel sizes', subtitle='PTA neurons only, |size|<20 bp') +
    ggtheme + theme(aspect.ratio=2/3)



# B. Fraction of 2bp dels among all indels
muts.new.ol.phenotype <- copy(muts)[, phenotype := fct_relevel(factor(ifelse(celltype=='oligo', 'Control oligo', ifelse(celltype=='neuron' & phenotype=='Control', 'Control neuron', as.character(phenotype)))), c('Control neuron', 'Control oligo', 'ALS', 'FTD', 'Alzheimers'))]
fig3b <- ggplot(muts.new.ol.phenotype[age > 50 & muttype == 'indel' & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & amp == 'PTA', .(frac.2bp=mean(nchar(refnt)-nchar(altnt)==2)), by=.(sample, phenotype)], aes(x=fct_drop(phenotype), y=frac.2bp, fill=fct_drop(phenotype))) +
    geom_boxplot(outliers=FALSE) +
    geom_jitter(size=1/4, width=0.2) +
    ylab('#2 bp deletions / all sIndels') + xlab(element_blank()) +
    labs(title='Indel sizes', subtitle='PTA neurons only, age>50, no MAPD/indel outliers') +
    ggtheme + theme(aspect.ratio=5/4) +
    expand_limits(y=0:1) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    geom_signif(comparisons=list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')), margin_top=-0.35, textsize=2.5, step_increase=1/10, size=1/4)



# C. De novo spectra and COSMIC ID4 for reference
id.sigs <- fread('mutsigs/sigprofilerextractor/all/ID83_corrected/Most_Stab_Sigs/signatures.csv')
# convert COSMIC format to SCAN2 matrix format
id.sigs.mat <- as.matrix(id.sigs, rownames=1)[names(table(id83(c()))),]
# Add COSMIC ID4
id.sigs.mat <- cbind(id.sigs.mat, as.matrix(fread('v3.3_ID4_PROFILE.txt'), rownames=1))
colnames(id.sigs.mat) <- c('ID-A', 'ID-B', 'COSMIC ID4')
fig3c <- wrap_elements(panel=~plot.id83(id.sigs.mat, show.detailed.types=F, uniform.y.axis=F, las=2), clip=FALSe)



# D. Boxplots of ID-A/B burdens vs. phenotype
fig3d <- ggplot(resids[burden.source=='denovo' & burden.type=='ID83_corrected' & age > 50 & sample %in% no.mapd.indel.outliers & tissue.origin != 'DG' & celltype == 'neuron' & amp == 'PTA'], aes(x=phenotype, y=pmax(1/4, (1+burden)/(1+expected.burden)), fill=phenotype)) +
    geom_boxplot(outlier.shape=NA, col=1, linewidth=0.3) +
    geom_jitter(col=1, size=1/2) +
    facet_wrap(~burden.subtype, ncol=3) +
    ggtheme +
    xlab('') + ylab('Signature exposure (obs / exp)') +
    guides(fill=guide_legend(title='Phenotype')) +
    geom_hline(yintercept=1, linewidth=0.5) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), aspect.ratio=2.0) +
    geom_signif(comparisons=list(c('Control', 'ALS'), c('Control', 'FTD'), c('Control', 'Alzheimers')),
        step_increase=1/10, margin_top=-0.1, textsize=2.75, size=0.15) +
    labs(title='De novo ID83 burdens by phenotype', subtitle='MAPD/indel outliers removed, age>50, floor set=0.25') +
    scale_y_continuous(trans='log2', n.breaks=8, label=function(x) round(x,2))


# E. Fraction of pink peak positive cells with odds ratio test
class.tab2 <- meta[amp=='PTA' & sample %in% no.mapd.indel.outliers,.(`ID-A high`=sum(delsig=='high'|delsig=='normal'), `ID-A normal`=sum(delsig=='negative'), total=nrow(.SD), rate=mean(delsig!='negative')),by=.(phenotype=factor(ifelse(phenotype=='Control',ifelse(celltype=='oligo', 'Control oligo', 'Control neuron'), as.character(phenotype)), levels=c('Control neuron', 'Control oligo', 'ALS', 'FTD', 'Alzheimers')))]

fishertests2 <- rbindlist(lapply(list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')), function(comp) {
    m <- as.matrix(class.tab2[phenotype %in% comp, .(`ID-A normal`, `ID-A high`)])  # normal first
    ft <- fisher.test(m)
    ci <- ft$conf.int
    data.frame(pheno1=comp[1], pheno2=comp[2],
        pheno1.n=sum(m[1,]), pheno1.positive=m[1,'ID-A high'],
        pheno2.n=sum(m[2,]), pheno2.positive=m[2,'ID-A high'],
        p.value=ft$p.value, odds.ratio=ft$estimate, ci.lower=ci[1], ci.upper=ci[2], ci.level=attr(ci, 'conf.level'))        
}))
fishertests2[, pheno1.frac := pheno1.positive/pheno1.n]
fishertests2[, pheno2.frac := pheno2.positive/pheno2.n]

# R's fisher.test() will not allow adding 0.5, so compute it manually
# So the OR here will not exactly match the CI computed by fisher.test
fishertests2[, odds.ratio.haldane.anscombe := ((pheno2.positive+0.5)/(pheno2.n-pheno2.positive+0.5))/((pheno1.positive+0.5)/(pheno1.n-pheno1.positive+0.5))]

fishertests2[, pheno2 := factor(pheno2, levels=c('Control oligo', 'ALS', 'FTD', 'Alzheimers'))]

fig3e <- ggplot(class.tab2, aes(x=phenotype, weight=rate, fill=phenotype, label=paste(`ID-A high`, total, sep='/'))) +
    geom_bar() +
    geom_text(aes(y=rate), nudge_y=0.035, size=3) +
    ggtheme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), aspect.ratio=5/4) +
    labs(title='Rates of ID-A high neurons', subtitle='PTA only, all ages, no MAPD/indel outliers') +
    xlab(element_blank()) + ylab('Fraction ID-A high cells') +
    expand_limits(y=c(0,1)) +
    geom_signif(aes(y=rate), comparisons=list(c('Control neuron', 'Control oligo'), c('Control neuron', 'ALS'), c('Control neuron', 'FTD'), c('Control neuron', 'Alzheimers')),
        step_increase=1/12, margin_top=-0.21, textsize=2.75, size=0.15, annotation=sprintf('%0.2g', fishertests2$p.value), position=position_nudge(y=-0.05))


# F. Spectra by cell type/phenotype and pink peak class
muts2 <- muts[amp == 'PTA' & muttype == 'indel' & sample %in% no.mapd.indel.outliers]
muts2.split <- split(muts2, paste(ifelse(muts2$delsig=='negative', 'negative', 'positive'), muts2$phenotype, muts2$celltype))
muts2.split <- muts2.split[c(1:3, 5, 4, 6:9)]  # hack to reorder the control oligo panel
fig3f <- plot.id83(sapply(muts2.split, function(df) table(id83(df$mutsig))), max.nrow=5, uniform.y.axis=F)




#
# 28. Gel single strand smear intensity quantitation (ImageJ) by phenotype
#
dmeta <- fread('metadata/immutable_donor_metadata.csv')
dmeta[, phenotype := factor(fct_recode(phenotype, AD='Alzheimers'), c('Control', 'FTD', 'AD'))]
ssi <- dmeta[fread('single_strand_intensity_gel_values.csv'), , on=.(donor)]
ctrl.mean <- ssi[phenotype == 'Control', mean(intensity)]
ssi[, intensity := intensity / ctrl.mean]

ggplot(ssi, aes(x=phenotype, y=intensity, fill=phenotype)) +
    geom_boxplot(outliers=F) + geom_jitter(width=0.2) +
    ggtheme + theme(aspect.ratio=3/2) +
    xlab(element_blank()) + ylab('Single-strand genomic breaks\n(relative intensity, AU)') +
    geom_signif(comparisons=list(c('Control', 'FTD'), c('Control', 'AD')), step_increase=1/12) +
    scale_y_continuous(expand=expansion(add=c(1, 1.)))
dev.print(dev=pdf, file='single_strand_break_intensity_compare_phenotypes.pdf')

# alternative: do a barplot of the mean with jittered points
#geom_bar(stat='summary', fun=mean)


#
# 29. SSB intensity vs. ID4 level
#
ggplot(resids[burden.type=='ID83_corrected' & burden.subtype=='ID83A' & burden.source=='denovo' & amp=='PTA' & celltype=='neuron' & sample %in% no.mapd.indel.outliers & phenotype != 'Control', .(mean.IDA=mean(obs.exp, na.rm=TRUE)), by=.(donor)][ssi, , on=.(donor)], aes(x=mean.IDA, y=intensity, col=fct_drop(phenotype))) +
    geom_point() + geom_smooth(method='lm', se=F) +
    ggtheme + theme(aspect.ratio=1) +
    ggpmisc::stat_fit_glance(method='lm', aes(label=sprintf('P = %0.2g, R^2 = %0.2g', after_stat(p.value), after_stat(r.squared))), size=3.5, label.x='right', label.y='bottom') +
    xlab('ID-A exposure') +
    ylab('Single-strand genomic breaks\n(relative intensity, AU)')
dev.print(dev=pdf, file='single_strand_break_intensity_vs_ID-A_exposure.pdf')



#
# 30. SSB intensity vs. age
#
ssi[, ida.high := TRUE ]
# These are taken from Zinan's plots, which were defined as "majority ID-A-high classified"
# It would be nice to compute this here to record the method.
# phenotype=='Control' : these were determined to be ID-A low by the method, despite looking
# like they're just declared to be that way.
ssi[donor %in% c('HBTRCS08631', 'MADRC1456', 'MADRC2207') | phenotype == 'Control', ida.high := FALSE]
ggplot(ssi, aes(x=age, y=intensity, col=fct_drop(phenotype), shape=ida.high)) +
    geom_point(size=2) +
    ggtheme + theme(aspect.ratio=1) +
    xlab('Age') +
    ylab('Single-strand genomic breaks\n(relative intensity, AU)') +
    scale_shape_manual(values=c(16,17))
dev.print(dev=pdf, file='single_strand_break_intensity_vs_age.pdf')




##############################################################################
# Supplemental Tables
##############################################################################

#
# Supplemental Table 1 - counts of cells and individuals
#
table1 <- meta[,.(age=age[1], sex=sex[1], phenotype=phenotype[1]), by=.(donor)]
fwrite(table1, file='table1.csv', quote=FALSE)

#
# Supplemental Table 2 - numbers of cells passing QC criteria
#
table2 <- meta[amp != 'bulk',.(`Mean age`=round(mean(age),1), Individuals=length(unique(donor)), `N male cells`=sum(sex=='M'), Cells=length(unique(sample)), `TDP43- cells`=sum(tdp43 =='Depleted'), `PFC cells`=sum(tissue.origin=='PFC'), `BA6 cells`=sum(tissue.origin=='BA6'), `DG cells`=sum(tissue.origin=='DG')),by=.(amp, phenotype, celltype, QC=ifelse(sample %in% no.mapd.indel.outliers, 'PASS', 'FAIL'))]
fwrite(table2, file='table2.csv', quote=FALSE)

#
# Supplemental Table 3 - detailed list of single cells with metadata
#

table3 <- fread('metadata/sample_metadata.csv')
table3 <- table3[amp != 'bulk']   # remove bulk IDs
# rename to match language in text/figures
table3[, delsig := ifelse(delsig=='negative', 'ID-A normal', 'ID-A high')]
table3[, phenotype := ifelse(phenotype == 'Alzheimers', 'AD', phenotype)]
# assign a numerical ID. this is only done to help reduce the size of
# the mutation call table, which annotates each call with the sample it
# was called in. recording the full sample name simply wastes too much space.
table3[, sample.id := 1:nrow(table3)]
fwrite(table3, file='table3.csv')


#
# Supplemental Table 4 - unfiltered mutation call set. Filtered calls are marked
# N.B requires table3 to reassign sample names to numerical IDs
#
table4 <- fread('tables/all___UNFILTERED_mut___any.csv')
table4 <- table4[chr %in% 1:22]    # remove X, Y calls
table4 <- table4[rescue == FALSE]  # no mutational signature rescue

# The file is very large, perhaps surpassing file size limits. Remove
# columns that aren't absolutely necessary.
table4[, id := NULL]
table4[, nearest := NULL]
table4[, pass := NULL]    # no longer necessary since rescued calls removed
table4[, rescue := NULL]  # no longer necessary since rescued calls removed

# only retain the "final.filter" column to save space. These 3 separate
# filters indicate why the final.filter is set.
table4[, rec.filter := NULL]
table4[, cluster.filter := NULL]
table4[, lineage.filter := NULL]

table4[, subject := NULL]

# translate sample name strings to numerical sample IDs and remove
# the (long) sample names to save space
table4 <- table3[,.(sample, sample.id)][table4,,on=.(sample)]
table4[, sample := NULL]

# rename the 'final.filter' column something more intelligible
setnames(table4, 'final.filter', 'fails.post.SCAN2.filters')

fwrite(table4, file='table4.csv') #, logical01=TRUE)


#
# Supplemental Table 5 - matrix of mutation/signature burdens
#

# Remove direct COSMIC fits and ID83 prior to correction for SCAN2
# sensitivity (ID83_corrected is what we want).  Rename columns
# to be more intelligible.
table5 <- resids[burden.source != 'cosmic' & burden.type != 'ID83', .(sample, type=burden.subtype, `mutation.count`=num.muts, `genome.wide.burden`=burden, expected.burden, obs.exp=(1+burden)/(1+expected.burden))][order(sample, type)]
fwrite(table5, file='table5.csv')
