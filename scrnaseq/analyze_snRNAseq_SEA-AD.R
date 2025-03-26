#### Analyze SEA-AD 10X snRNA-seq data ####

options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(schard)

# Load h5ad as Seurat by cell type (24 total) ####
out_dir <- "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass"
file_list <- list.files(out_dir, "h5ad")
for (file in file_list){
    file = file_list[i]
    subclass <- unlist(str_split(file, "\\."))[4]
    print(subclass)
    # Convert h5ad to seurat object
    seurat_obj <- schard::h5ad2seurat(paste0(out_dir, "/", file))
    # Replace counts with UMI count matrix
    umi <- readMM(str_replace(paste0(out_dir, "/", file),"h5ad","UMIs.mtx"))
    umi_new <- as(t(umi), "CsparseMatrix")
    rownames(umi_new) <- rownames(seurat_obj[["RNA"]]$counts)
    colnames(umi_new) <- colnames(seurat_obj[["RNA"]]$counts)
    seurat_obj[["RNA"]]$counts <- umi_new
    # Save seurat object
    saveRDS(seurat_obj, str_replace(paste0(out_dir, "/", file),"h5ad","rds"))
}
# Load and convert L2/3 IT by 10 supertypes ####
out_dir <- "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass"
file_list <- list.files(out_dir, "h5ad")
file_list <- file_list[grepl("L2_3_IT_[0-9]", file_list)]
for (file in file_list){
    supertype <- unlist(str_split(file, "\\."))[4]
    print(supertype)
    # Convert h5ad to seurat object
    seurat_obj <- schard::h5ad2seurat(paste0(out_dir, "/", file))
    # Replace counts with UMI count matrix
    umi <- readMM(str_replace(paste0(out_dir, "/", file),"h5ad","UMIs.mtx"))
    umi_new <- as(t(umi), "CsparseMatrix")
    rownames(umi_new) <- rownames(seurat_obj[["RNA"]]$counts)
    colnames(umi_new) <- colnames(seurat_obj[["RNA"]]$counts)
    seurat_obj[["RNA"]]$counts <- umi_new
    # Save seurat object
    saveRDS(seurat_obj, str_replace(paste0(out_dir, "/", file),"h5ad","rds"))
}
# Load and convert total L2/3 IT after filtering for DEG ####
out_dir <- "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass"
file <- "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.L2_3_IT.DEG_AMPAD.h5ad"
# Convert h5ad to seurat object
seurat_obj <- schard::h5ad2seurat(paste0(out_dir, "/", file))
# Replace counts with UMI count matrix
umi <- readMM(str_replace(paste0(out_dir, "/", file),"h5ad","UMIs.mtx"))
umi_new <- as(t(umi), "CsparseMatrix")
rownames(umi_new) <- rownames(seurat_obj[["RNA"]]$counts)
colnames(umi_new) <- colnames(seurat_obj[["RNA"]]$counts)
seurat_obj[["RNA"]]$counts <- umi_new
# Save seurat object
saveRDS(seurat_obj, str_replace(paste0(out_dir, "/", file),"h5ad","rds"))


# SEA-AD DLPFC DEG -------------------------------------------------------------------
library(nebula)
out_dir <- "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass"
file_list <- list.files(out_dir, "h5ad")
for (file in file_list){
    supertype <- unlist(str_split(file, "\\."))[4]
    print(supertype)
    seurat_obj <- readRDS(str_replace(paste0(out_dir, "/", file),"h5ad","rds"))
    seurat_obj$nCount_RNA <- seurat_obj$Number_of_UMIs
    # Create variables used in model.matrix
    seurat_obj$Race <- case_when(
        seurat_obj$Race_choice_Other == "Checked" ~ "Other",
        seurat_obj$Race_choice_White == "Checked" ~ "White",
        seurat_obj$Race_choice_Black_African_American == "Checked" ~ "Black_African_American",
        seurat_obj$Race_choice_Asian == "Checked" ~ "Asian",
        seurat_obj$Race_choice_American_Indian_Alaska_Native == "Checked" ~ "American_Indian_Alaska_Native",
        seurat_obj$Race_choice_Native_Hawaiian_or_Pacific_Islander == "Checked" ~ "Native_Hawaiian_or_Pacific_Islander",
        seurat_obj$Race_choice_Unknown_or_unreported == "Checked" ~ "Unknown_or_unreported"
        )
    meta <- read.csv("AD_meta_Table_S4_mCA_sample_list_individual_metadata.csv",header=T)
    meta <- meta[meta$Individual.ID!="",1:7]
    control_donors <- meta$Individual.ID[meta$Diagnosis=="Healthy control"]
    ad_donors <- meta$Individual.ID[meta$Diagnosis=="Alzheimer's disease"]
    seurat_obj <- subset(seurat_obj, Donor_ID%in%meta$Individual.ID)
    seurat_obj$DEG_group <- "Unknown"
    seurat_obj$DEG_group[seurat_obj$Donor_ID%in%control_donors] <- "Control"
    seurat_obj$DEG_group[seurat_obj$Donor_ID%in%ad_donors] <- "AD"
    # Filter out Severely Affected donors and unknown DEG_group
    seurat_dge <- subset(seurat_obj, Severely_Affected_Donor=="N" & DEG_group!="Unknown")
    seurat_dge$DEG_group <- factor(seurat_dge$DEG_group, levels=c("Control","AD"))
    # Convert seurat to nebula 
    seuratdata <- scToNeb(obj = seurat_dge, assay = "RNA", id = "Donor_ID", pred = c("DEG_group", "Sex", "Age_at_Death", "Race", "method", "Genes_detected"), offset="nCount_RNA")
    # Build model matrix
    df = model.matrix(~ DEG_group, data=seuratdata$pred)
    # Add pseudocount
    data_g = group_cell(count=seuratdata$count, id=seuratdata$id, pred=df, offset=seuratdata$offset)
    re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset, cpc=0.005)
    saveRDS(re, str_replace(paste0(out_dir, "/", file),"h5ad","DEG_Dementia.nebula.rds"))
}

# Heatmap of DEG results ####
target_genes <- c("RRM1", "RRM2B", "RNASEH2A", "RNASEH2B", "RNASEH2C")
files <- list.files("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/", "nebula")
celltypes <- str_remove_all(str_remove_all(files, "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean."), ".DEG_AMPAD.nebula.rds")
df <- c()
for (i in 1:length(files)){
    f <- files[i]
    celltype <- celltypes[i]
    tmp <- readRDS(paste0("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/",f))
    tmp <- tmp$summary
    tmp$p_val_adj <- p.adjust(tmp$p_DEG_groupAD, method = "BH")
    tmp <- tmp[tmp$gene %in% target_genes,]
    tmp$avg_log2FC <- tmp$logFC_DEG_groupAD
    tmp <- tmp[,c("gene","avg_log2FC","p_val_adj")]
    tmp$celltype <- celltype
    df <- rbind(df, tmp)
}
df.plot <- df
# Plot heatmap
target_genes <- c("RRM1", "RRM2B", "RNASEH2A", "RNASEH2B", "RNASEH2C")
celltypes <- sort(unique(df.plot$celltype))
ordered_celltypes <- celltypes[c(4:24,28:34)]
df.plot <- df.plot[df.plot$gene%in%target_genes&df.plot$celltype%in%ordered_celltypes,]
df.plot$mylabel <- ifelse(df.plot$p_val_adj<0.05, "*", "")
df.plot$gene <- factor(df.plot$gene, levels = rev(target_genes))
df.plot$celltype <- factor(df.plot$celltype, levels = ordered_celltypes)
p <- ggplot(df.plot, aes(x=celltype, y=gene, fill=avg_log2FC)) +
    geom_tile() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    geom_text(aes(label=mylabel)) +
    labs(title="SEA-AD DLPFC AD vs. Control (Nebula)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
p

