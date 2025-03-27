# Load necessary libraries
library(scDblFinder)        # Doublet detection using statistical modeling
library(DropletQC)          # Estimation of nuclear fraction from BAM and barcode
library(SingleCellExperiment) # Required for interoperability with scDblFinder and DropletQC
library(Seurat)             # Core framework for single-cell RNA-seq analysis (normalization, PCA, clustering, etc.)
library(harmony)            # Batch effect correction using Harmony
library(ggplot2)            # Visualization of UMAPs, violin plots, etc.
library(ggpubr)             # Enhancements for ggplot2 (e.g., easy publication-quality plots)
library(Matrix)             # Reading and manipulating sparse matrices (e.g., CellRanger outputs)
library(dplyr)              # Data wrangling, especially for metadata manipulation
library(BiocParallel)       # Parallel computing support, used in scDblFinder and DropletQC
library(patchwork)          # Combining multiple ggplot2 plots (e.g., UMAPs at different resolutions)

# Step 1: Set up parallel backend using SLURM environment variable
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
param <- MulticoreParam(workers = num_cores)
register(param)

print("This is the number of cores")
print(num_cores)

# Step 2: Define function to create Seurat object from sparse matrix
create_seurat_object_from_matrix <- function(group_dir, sample_name, original_or_filtered) {
  base_name <- paste0(sample_name, "_", original_or_filtered)
  mtx_file <- file.path(group_dir, sample_name, paste0(base_name, ".mtx"))
  barcodes_file <- file.path(group_dir, sample_name, paste0(base_name, "_barcodes.tsv"))
  genes_file <- file.path(group_dir, sample_name, paste0(base_name, "_genes.tsv"))
  counts <- readMM(mtx_file)
  barcodes <- readLines(barcodes_file)
  genes <- readLines(genes_file)
  rownames(counts) <- barcodes
  colnames(counts) <- genes
  counts <- t(counts)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  return(seurat_obj)
}

# Step 3: Add metadata information to Seurat object
add_metadata_our_group <- function(seurat_obj, sample_name, metadata_df) {
  metadata <- metadata_df[metadata_df$our_sample_name == sample_name, ]
  seurat_obj$case_status <- metadata$case_status
  seurat_obj$case_ID <- metadata$case_ID
  seurat_obj$brain_region <- metadata$brain_region
  seurat_obj$TDP43_status <- metadata$TDP43_status
  seurat_obj$condition_1 <- metadata$condition_1
  seurat_obj$condition_2 <- metadata$condition_2
  return(seurat_obj)
}

# Step 4: Read metadata and create Seurat objects for all samples
our_group_meta_df <- read.delim("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/data/our_group/sample_name.tsv", stringsAsFactors = FALSE)
our_group_dir <- "/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/11_cellbender_QC_sparse_matrix/our_group"
original_or_filtered <- "filtered" # Cellbender filtered matrix
sample_names <- list.dirs(our_group_dir, full.names = FALSE, recursive = FALSE)

seurat_objects <- list()
for (sample_name in sample_names) {
  seurat_obj <- create_seurat_object_from_matrix(our_group_dir, sample_name, original_or_filtered)
  seurat_obj <- add_metadata_our_group(seurat_obj, sample_name, our_group_meta_df)
  object_name <- paste0(sample_name, "_", original_or_filtered)
  assign(object_name, seurat_obj)
  seurat_objects[[object_name]] <- seurat_obj
}

# Step 5: Doublet detection and nuclear fraction estimation per sample
for (name in names(seurat_objects)) {
  print(paste("Processing:", name))
  seurat_object <- seurat_objects[[name]]
  sce <- as.SingleCellExperiment(seurat_object)
  sce <- scDblFinder(sce, BPPARAM = param)
  seurat_object@meta.data$doublet_class <- colData(sce)$scDblFinder.class
  sample_name <- strsplit(name, "_")[[1]][1]
  test_bam_path <- paste0('/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/01_cell_ranger_output_our_group/', sample_name, '/outs/possorted_genome_bam.bam')
  test_barcs <- rownames(seurat_object@meta.data)
  nf <- nuclear_fraction_tags(bam = test_bam_path, barcodes = test_barcs, tiles = 1, cores = num_cores - 1, verbose = TRUE)
  seurat_object$nuclear_fraction <- nf[match(rownames(seurat_object@meta.data), rownames(nf)), "nuclear_fraction"]
  seurat_objects[[name]] <- seurat_object
}

# Step 6: Merge all sample Seurat objects into one
merged_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)
saveRDS(merged_seurat_object, file = "/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_seurat_object.rds")

# Step 7: Visualize nuclear fraction distribution by UMI count range (before filtering)
seurat_object <- readRDS("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_seurat_object.rds")
seurat_object@meta.data <- seurat_object@meta.data %>%
  mutate(UMI_range = cut(nCount_RNA, breaks = c(seq(100, 2000, by = 100), Inf), right = FALSE, labels = paste(seq(100, 2000, by = 100), seq(200, 2100, by = 100), sep = "-"))) %>%
  mutate(UMI_range = factor(UMI_range, levels = paste(seq(100, 2000, by = 100), seq(200, 2100, by = 100), sep = "-")))
nuclear_fraction_plot <- ggplot(seurat_object@meta.data, aes(x = UMI_range, y = nuclear_fraction, fill = UMI_range)) +
  geom_violin() +
  ggtitle("Distribution of Nuclear Fraction by UMI Count Range") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.text.x = element_text(size = 20, angle = 40, hjust = 1), axis.text.y = element_text(size = 12)) +
  xlab("UMI Count Range") + ylab("Nuclear Fraction")
print(nuclear_fraction_plot)

# Step 8: Filter by sample, doublets, and nuclear fraction threshold
filtered_seurat_object <- subset(seurat_object, subset = orig.ident != "1700BA6-normal") # lower_qulaity
doublet_filtered_seurat_object <- subset(filtered_seurat_object, subset = doublet_class != "doublet")
NF_cut_off_val = 0.25
final_filtered_seurat_object <- subset(doublet_filtered_seurat_object , subset =  nuclear_fraction>=NF_cut_off_val)

# Step 9: Visualize nuclear fraction after filtering
final_filtered_seurat_object@meta.data <- final_filtered_seurat_object@meta.data %>%
  mutate(UMI_range = cut(nCount_RNA, breaks = c(seq(100, 2000, by = 100), Inf), right = FALSE, labels = paste(seq(100, 2000, by = 100), seq(200, 2100, by = 100), sep = "-"))) %>%
  mutate(UMI_range = factor(UMI_range, levels = paste(seq(100, 2000, by = 100), seq(200, 2100, by = 100), sep = "-")))
nuclear_fraction_plot <- ggplot(final_filtered_seurat_object@meta.data, aes(x = UMI_range, y = nuclear_fraction, fill = UMI_range)) +
  geom_violin(drop = FALSE) +
  ggtitle("Distribution of Nuclear Fraction by UMI Count Range") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.text.x = element_text(size = 20, angle = 40, hjust = 1), axis.text.y = element_text(size = 12)) +
  xlab("UMI Count Range") + ylab("Nuclear Fraction")
print(nuclear_fraction_plot)
saveRDS(final_filtered_seurat_object, file = "/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_final_QC_seurat_object.rds")

# Step 10: SCTransform normalization and PCA
final_filtered_seurat_object = readRDS("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_final_QC_seurat_object.rds")
final_filtered_seurat_object <- PercentageFeatureSet(final_filtered_seurat_object, pattern = "^MT-", col.name = "percent.mt")
final_filtered_seurat_object <- SCTransform(final_filtered_seurat_object, vars.to.regress = "percent.mt", verbose = TRUE, vst.flavor = "v2")
final_filtered_seurat_object <- RunPCA(final_filtered_seurat_object, npcs = 30, verbose = TRUE)
saveRDS(final_filtered_seurat_object, file = "/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_SCTtransform.rds")

# Step 11: Determine optimal number of PCs
final_filtered_seurat_object <- readRDS("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_SCTtransform.rds")
pct <- final_filtered_seurat_object[["pca"]]@stdev / sum(final_filtered_seurat_object[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# Step 12: Batch effect correction using Harmony
final_filtered_seurat_object <- RunHarmony(final_filtered_seurat_object, group.by.vars = "orig.ident", dims.use = 1:17)
saveRDS(final_filtered_seurat_object, file = "/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_harnomy_with_17_PCs.rds")

# Step 13: Clustering and UMAP visualization
final_filtered_seurat_object <- readRDS("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/12_seurat_object/our_sample/gpu_merged_after_harnomy_with_17_PCs.rds")
final_filtered_seurat_object <- FindNeighbors(final_filtered_seurat_object, reduction = "harmony", dims = 1:17)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
umap_plots <- list()
for (res in resolutions) {
  final_filtered_seurat_object <- FindClusters(final_filtered_seurat_object, resolution = res, algorithm = 4, method = "igraph", verbose = TRUE)
  final_filtered_seurat_object <- RunUMAP(final_filtered_seurat_object, reduction = "harmony", dims = 1:17)
  umap_plot <- DimPlot(final_filtered_seurat_object, reduction = "umap", group.by = "seurat_clusters") + ggtitle(paste("Resolution:", res))
  umap_plots[[as.character(res)]] <- umap_plot
}
combined_umap_plot <- wrap_plots(umap_plots, ncol = 2) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + theme(legend.position = "bottom")
ggsave("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/code/20_new_downstream_analysis/07_clustering/our_group_leiden/combined_umap_plots.pdf", plot = combined_umap_plot, width = 20, height = 28)

# Step 14: Visualize clustering consistency using clustree
clustree_plot <- clustree(final_filtered_seurat_object)
ggsave("/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/code/20_new_downstream_analysis/07_clustering/our_group_leiden/clustree_plot.pdf", plot = clustree_plot, width = 28, height = 16)
