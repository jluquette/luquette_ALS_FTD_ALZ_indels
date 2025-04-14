library(Seurat)
library(dplyr)
library(nebula)
library(BiocParallel)

# Get the number of available cores from SLURM environment variables
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# Set up parallel backend
param <- MulticoreParam(workers = num_cores)
register(param)
print("this is number of cores")
print(num_cores)

# Load Seurat object of Manolis data's their clean up version
##For each differential expression cluster, only genes present in at least 10% of cells were retained for analysis
cat("Loading Seurat object\n")
our_group_so<- readRDS("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/Final_three/seurat_object/gpu_merged_after_harmony_with_17_PCs.cluster_ready.marker_ready.annotated.07152024.rds")
DefaultAssay(our_group_so) <- "RNA"
cat("Below is the seurat_object\n")
print(our_group_so)
head(our_group_so@meta.data)


## Case 1 :  offset = "nCount_RNA", term = ~ case_status * nuclear_fraction )
# Loop through each unique cell type and perform analysis for both FTD vs. Normal(Ctrl) and ALS vs. Normal(Ctrl)
our_group_so@meta.data$celltype_detailed <- as.character(our_group_so@meta.data$celltype_detailed)
celltypes <- unique(our_group_so@meta.data$celltype_detailed)

# Initialize lists to store results
results_ln_ftd_list <- list()
results_ln_als_list <- list()

for (celltype in celltypes) {
    cat("Processing cell type:", celltype, "\n")

    # Subset data for the specific cell type
    #subset_obj <- subset(our_group_so, subset = celltype_detailed == celltype)

    # Subset data for the specific cell type
    celltype_indices <- which(our_group_so@meta.data$celltype_detailed == celltype)
    subset_obj <- our_group_so[, celltype_indices]
    DefaultAssay(subset_obj) <- "RNA"
    cat("Subset based on cell type : done", "\n")

    ### FTD vs. Normal(Ctrl) ###
    specific_case_status <- c("FTD","Normal(Ctrl)")
    specific_brain_region <- "PFC"
    subset_obj_ftd <- subset(subset_obj,subset = case_status %in% specific_case_status & brain_region == specific_brain_region)
    #subset_obj_ftd <- subset(subset_obj, subset = case_status %in% specific_case_status)
    subset_obj_ftd@meta.data$case_status <- factor(subset_obj_ftd@meta.data$case_status, levels = specific_case_status)
    cat("Subset based on case_status : done", "\n")

    seuratdata_ftd <- scToNeb(obj = subset_obj_ftd, assay = "RNA", id = "case_ID", pred = c("case_status", "nuclear_fraction"), offset = "nCount_RNA")
    df_ftd <- model.matrix(~ case_status * nuclear_fraction, data = seuratdata_ftd$pred)
    data_g_ftd <- group_cell(
        count = seuratdata_ftd$count,
        id = seuratdata_ftd$id,
        pred = df_ftd,
        offset = seuratdata_ftd$offset
    )

    if (is.null(data_g_ftd)) {
        re_ln_ftd <- nebula(seuratdata_ftd$count, seuratdata_ftd$id, pred = df_ftd, offset = seuratdata_ftd$offset, method = 'LN', ncore = num_cores - 1)
    } else {
        re_ln_ftd <- nebula(data_g_ftd$count, data_g_ftd$id, pred = df_ftd, offset = data_g_ftd$offset, method = 'LN', ncore = num_cores - 1)
    }

    # Store FTD vs. Normal(Ctrl) results in the respective lists
    results_ln_ftd_list[[celltype]] <- re_ln_ftd
    ### ALS vs. Normal(Ctrl) ###
    specific_case_status <- c("ALS", "Normal(Ctrl)")
    specific_brain_region <- "BA6"
    subset_obj_als <- subset(subset_obj,subset = case_status %in% specific_case_status & brain_region == specific_brain_region)
    #subset_obj_als <- subset(subset_obj, subset = case_status %in% specific_case_status)
    subset_obj_als@meta.data$case_status <- factor(subset_obj_als@meta.data$case_status, levels = specific_case_status)

    seuratdata_als <- scToNeb(obj = subset_obj_als, assay = "RNA", id = "case_ID", pred = c("case_status", "nuclear_fraction"), offset = "nCount_RNA")
    df_als <- model.matrix(~ case_status * nuclear_fraction, data = seuratdata_als$pred)
    data_g_als <- group_cell(
        count = seuratdata_als$count,
        id = seuratdata_als$id,
        pred = df_als,
        offset = seuratdata_als$offset
    )

    if (is.null(data_g_als)) {
        re_ln_als <- nebula(seuratdata_als$count, seuratdata_als$id, pred = df_als, offset = seuratdata_als$offset, method = 'LN', ncore = num_cores - 1)
    } else {
        re_ln_als <- nebula(data_g_als$count, data_g_als$id, pred = df_als, offset = data_g_als$offset, method = 'LN', ncore = num_cores - 1)
    }

    # Store ALS vs. Normal(Ctrl) results in the respective lists
    results_ln_als_list[[celltype]] <- re_ln_als
}

# Save the results to .rds files for FTD vs. Normal(Ctrl)
saveRDS(results_ln_ftd_list, file = "./03_results_our_group/our_results_ln_list_FTD_vs_Normal_IT.rds")
cat("LN results for FTD vs. Normal(Ctrl) have been saved to ./03_results_our_group/our_results_ln_list_FTD_vs_Normal_IT.rds\n")

# Save the results to .rds files for ALS vs. Normal(Ctrl)
saveRDS(results_ln_als_list, file = "./03_results_our_group/our_results_ln_list_ALS_vs_Normal_IT.rds")
cat("LN results for ALS vs. Normal(Ctrl) have been saved to ./03_results_our_group/our_results_ln_list_ALS_vs_Normal_IT.rds\n")



# ### Case 2 :  no interaction term
# Loop through each unique cell type and perform analysis for both FTD vs. Normal(Ctrl) and ALS vs. Normal(Ctrl)
celltypes <- unique(our_group_so@meta.data$celltype_detailed)

# Initialize lists to store results
results_ln_ftd_list_no_IT <- list()
results_ln_als_list_no_IT <- list()

for (celltype in celltypes) {
    cat("Processing cell type:", celltype, "\n")

    # Subset data for the specific cell type
    #subset_obj <- subset(our_group_so, subset = celltype_detailed == celltype)
    # Subset data for the specific cell type
    celltype_indices <- which(our_group_so@meta.data$celltype_detailed == celltype)
    subset_obj <- our_group_so[, celltype_indices]
    DefaultAssay(subset_obj) <- "RNA"

    ### FTD vs. Normal(Ctrl) ###
    specific_case_status <- c("FTD", "Normal(Ctrl)")
    specific_brain_region <- "PFC"
    subset_obj_ftd <- subset(subset_obj,subset = case_status %in% specific_case_status & brain_region == specific_brain_region)
    #subset_obj_ftd <- subset(subset_obj, subset = case_status %in% specific_case_status)
    subset_obj_ftd@meta.data$case_status <- factor(subset_obj_ftd@meta.data$case_status, levels = specific_case_status)

    seuratdata_ftd <- scToNeb(obj = subset_obj_ftd, assay = "RNA", id = "case_ID", pred = c("case_status"), offset = "nCount_RNA")
    df_ftd <- model.matrix(~ case_status, data = seuratdata_ftd$pred)
    data_g_ftd <- group_cell(
        count = seuratdata_ftd$count,
        id = seuratdata_ftd$id,
        pred = df_ftd,
        offset = seuratdata_ftd$offset
    )

    if (is.null(data_g_ftd)) {
        re_ln_ftd <- nebula(seuratdata_ftd$count, seuratdata_ftd$id, pred = df_ftd, offset = seuratdata_ftd$offset, method = 'LN', ncore = num_cores - 1)
    } else {
        re_ln_ftd <- nebula(data_g_ftd$count, data_g_ftd$id, pred = df_ftd, offset = data_g_ftd$offset, method = 'LN', ncore = num_cores - 1)
    }

    # Store FTD vs. Normal(Ctrl) results in the respective lists
    results_ln_ftd_list_no_IT[[celltype]] <- re_ln_ftd

    ### ALS vs. Normal(Ctrl) ###
    specific_case_status <- c("ALS", "Normal(Ctrl)")
    specific_brain_region <- "BA6"
    subset_obj_als <- subset(subset_obj,subset = case_status %in% specific_case_status & brain_region == specific_brain_region)
    #subset_obj_als <- subset(subset_obj, subset = case_status %in% specific_case_status)
    subset_obj_als@meta.data$case_status <- factor(subset_obj_als@meta.data$case_status, levels = specific_case_status)

    seuratdata_als <- scToNeb(obj = subset_obj_als, assay = "RNA", id = "case_ID", pred = c("case_status"), offset = "nCount_RNA")
    df_als <- model.matrix(~ case_status, data = seuratdata_als$pred)
    data_g_als <- group_cell(
        count = seuratdata_als$count,
        id = seuratdata_als$id,
        pred = df_als,
        offset = seuratdata_als$offset
    )

    if (is.null(data_g_als)) {
        re_ln_als <- nebula(seuratdata_als$count, seuratdata_als$id, pred = df_als, offset = seuratdata_als$offset, method = 'LN', ncore = num_cores - 1)
    } else {
        re_ln_als <- nebula(data_g_als$count, data_g_als$id, pred = df_als, offset = data_g_als$offset, method = 'LN', ncore = num_cores - 1)
    }

    # Store ALS vs. Normal(Ctrl) results in the respective lists
    results_ln_als_list_no_IT[[celltype]] <- re_ln_als
}

# Save the results to .rds files for FTD vs. Normal(Ctrl)
saveRDS(results_ln_ftd_list_no_IT, file = "./03_results_our_group/our_results_ln_list_FTD_vs_Normal_no_IT.rds")
cat("LN results for FTD vs. Normal(Ctrl) have been saved to ./03_results_our_group/our_results_ln_list_FTD_vs_Normal_no_IT.rds\n")

# Save the results to .rds files for ALS vs. Normal(Ctrl)
saveRDS(results_ln_als_list_no_IT, file = "./03_results_our_group/our_results_ln_list_ALS_vs_Normal_no_IT.rds")
cat("LN results for ALS vs. Normal(Ctrl) have been saved to ./03_results_our_group/our_results_ln_list_ALS_vs_Normal.rds_no_IT\n")
