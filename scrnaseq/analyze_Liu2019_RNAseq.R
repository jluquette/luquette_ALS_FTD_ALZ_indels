#### Analyze Liu et al. 2019 bulk RNA-seq data ####

options(stringsAsFactors = FALSE)
library(reshape2)
library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)

# read count matrix
mat <- read.table("data/Liu2019_RNAseq/GSE126542_NeuronalNuclei_RNAseq_counts.txt", header=T, sep="\t", row.names=NULL)
names(mat)[1] <- "Symbol"

# Normalize by each sample's total count, i.e. RPM (do not need to normalize by gene length b/c comparing same genes)
tmp <- as.matrix(mat[,2:ncol(mat)])
tmp1 <- t(t(tmp)/(colSums(tmp)/1000000)) #RPM
mat_rpm <- data.frame(tmp1)
mat_rpm$Symbol <- mat$Symbol
mat_rpm <- mat_rpm[,c(ncol(mat_rpm),1:(ncol(mat_rpm)-1))]
# Extract inhibitory marker genes
target_genes <- c("GAD1", "GAD2")
df <- mat_rpm[mat_rpm$Symbol%in%target_genes, ]
df <- melt(df, id.vars = "Symbol", variable.name = "sample_name", value.name = "RPM")
df$sample_name <- as.character(df$sample_name)
df$sample_ID <- unlist(lapply(strsplit(df$sample_name, split="_"), `[[`, 1))
df$TDP43_status <- unlist(lapply(strsplit(df$sample_name, split="_"), `[[`, 2))
unique(df$TDP43_status) #"TN" "TP"
df.plot <- df
df.plot$TDP43_status <- ifelse(df.plot$TDP43_status=="TP", "positive", "negative")

# Plot
p <- ggplot(df.plot, aes(x=TDP43_status, y=RPM)) +
  facet_grid(~Symbol) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3) +
  stat_compare_means(method = "wilcox.test", paired = T, ref.group = "positive", method.args = list(alternative = "less"), label="p.format") +
  labs(title="one-tailed paired wilcoxon") +
  theme_classic()
p





