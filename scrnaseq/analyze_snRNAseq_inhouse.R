#### Analyze in-house ALS/FTD 10X snRNA-seq data ####

options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(cowplot)
library(clustree)
library(harmony)
library(reticulate)
library(ggrepel)
library(gridExtra)
library(ggpubr)
use_condaenv(condaenv = "py37", required = TRUE)
py_run_string('import leidenalg')


# 1. Clustering ####
seurat_integrated <- readRDS("gpu_merged_after_harmony_with_17_PCs.rds")
# Find clusters using harmony embeddings
resolutions = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4)
seurat_integrated <- FindNeighbors(object = seurat_integrated, reduction = "harmony")
seurat_integrated <- FindClusters(object = seurat_integrated, method="igraph", algorithm = 4, #use leiden algorithm
                                  resolution = resolutions)
seurat_integrated@meta.data %>% View()
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:17)
# Plot UMAP for each resolution
umap_plots <- list()
# Loop through each resolution and store the plot
for (res in resolutions) {
  # Generate UMAP plot and store in the list
  umap_plot <- DimPlot(seurat_integrated, reduction = "umap", 
                       group.by = paste0("SCT_snn_res.",res),
                       label = TRUE, label.size = 5, repel = F) + 
    ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
  umap_plots[[as.character(res)]] <- umap_plot
}
# Combine UMAP plots into one plot for comparison with adjusted layout
combined_umap_plot <- wrap_plots(umap_plots, ncol = 2) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust margins to increase space between plots
  theme(legend.position = "bottom")  # Move legend to the bottom
# Display the combined UMAP plot
print(combined_umap_plot)
# Check clustree to determine resolution
p <- clustree(seurat_integrated, prefix = "SCT_snn_res.")
p
# Check clustering by sample
res = 0.8
DimPlot(seurat_integrated, reduction = "umap", 
        group.by = paste0("SCT_snn_res.",res),
        split.by = "orig.ident", ncol = 4,
        label = TRUE, label.size = 5, repel = F) + 
  ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")


# 2. General cell type annotation ####
# use resolution = 0.8 
Idents(object = seurat_integrated) <- "SCT_snn_res.0.8"
# Known markers
markers <- list(
  neuronal = c("RBFOX3", "FOXP2", "SYNPR"),
  non_neuronal = c("GFAP", "OLIG2", "CSF1R", "CLDN5", "PTPRC"),
  excitatory_general = c("SLC17A6", "SLC17A7"),
  excitatory_upper = c("CUX2", "GLRA3", "LAMP5"),
  excitatory_middle = c("BHLHE22", "RBP4", "RORB"),
  excitatory_lower = c("NR4A2", "TLE4"),
  inhibitory_general = c("CCK", "GAD1", "GAD2", "SLC32A1"),
  inhibitory_subtypes = c("GAD1", "LAMP5", "VIP", "SST", "PVALB"),
  additional = c("TBR1", "NR4A2", "RBP4", "PAX6", "SOX6", "KCNK2"),
  Ch_IN = c("GAD1","UNC5B","RORA","TRPS1","NFIB"),
  L4_5_IT = c("RORB","GPR26","TMSB10"),
  L5_corticobulbar_tract = c("GRIK1","SLC24A3"), 
  L5_Betz_VEN = c("POU3F1"),
  L6b_CT = c("SULF1","SRGAP1","SEMA5A","SYT6","ZFHX3"),
  L6_IT_Car3 = c("NR4A2","SYT6","NTNG2","POSTN"),
  L6_IT = c("ITGB8","SNED1","OPRK1") 
)
# Annotate clusters with known markers and visualize
for (marker_group in names(markers)) {
  marker_plots <- lapply(markers[[marker_group]], function(marker) {
    p <- FeaturePlot(seurat_integrated,
                     reduction = "umap",
                     features = marker,
                     order = TRUE,
                     min.cutoff = 'q10',
                     label = TRUE) +
      ggtitle(marker)
    return(p)
  })
  combined_plot <- wrap_plots(marker_plots, ncol = 3) +
    plot_annotation(title = paste(marker_group, "markers")) & 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  print(combined_plot)
}

# Assign cell types ####
res <- "0.8"
Idents(object = seurat_integrated) <- paste0("SCT_snn_res.",res)
nclusters <- nlevels(seurat_integrated@active.ident)
##general cell type
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "16" = "Astrocyte",
                                  "13" = "Endothelial / Oligodentrocyte / Microglia / Blood",
                                  "1" = "Excitatory", "5" = "Excitatory", "2" = "Excitatory", "9" = "Excitatory", "8" = "Excitatory", "3" = "Excitatory", "4" = "Excitatory", "17" = "Excitatory", "18" = "Excitatory", "14" = "Excitatory", "11" = "Excitatory", "15" = "Excitatory",
                                  "6" = "Inhibitory", "12" = "Inhibitory", "19" = "Inhibitory", "7" = "Inhibitory", "10" = "Inhibitory")
seurat_integrated$celltype_general <- seurat_integrated@active.ident
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = F)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
seurat_integrated <- PrepSCTFindMarkers(object = seurat_integrated)
Idents(object = seurat_integrated) <- "SCT_snn_res.0.8"
markers_0.8 <- FindAllMarkers(seurat_integrated, assay = "SCT", only.pos = TRUE)
markers_0.8 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p <- DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend() + ggtitle("Res = 0.8")
p


# 3. Cell subtype annotation ####
# Subset to excitatory, inhibitory, and non-neuroanl, then find clusters
excitatory <- c(1, 5, 2, 9, 8, 3, 4, 17, 18, 14, 11, 15)
inhibitory <- c(6, 12, 19, 7, 10)
nonneuronal <- c(16, 13)
Idents(object = seurat_integrated) <- "SCT_snn_res.0.8"
seurat_ex <- subset(seurat_integrated, subset = SCT_snn_res.0.8 %in% excitatory)
seurat_in <- subset(seurat_integrated, subset = SCT_snn_res.0.8 %in% inhibitory)
seurat_non <- subset(seurat_integrated, subset = SCT_snn_res.0.8 %in% nonneuronal)

seurat_ex <- RunPCA(seurat_ex, npcs = 30, verbose = FALSE)
seurat_in <- RunPCA(seurat_in, npcs = 30, verbose = FALSE)
seurat_non <- RunPCA(seurat_non, npcs = 30, verbose = FALSE)
ElbowPlot(object = seurat_ex, ndims = 30)
ElbowPlot(object = seurat_non, ndims = 30)
# Determine percent of variation associated with each PC
pct <- seurat_ex[["pca"]]@stdev / sum(seurat_ex[["pca"]]@stdev) * 100
pct <- seurat_in[["pca"]]@stdev / sum(seurat_in[["pca"]]@stdev) * 100
pct <- seurat_non[["pca"]]@stdev / sum(seurat_non[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
# Run Harmony
seurat_ex <- RunHarmony(seurat_ex, group.by.vars = "orig.ident", dims.use = 1:25)
seurat_in <- RunHarmony(seurat_in, group.by.vars = "orig.ident", dims.use = 1:20)
seurat_non <- RunHarmony(seurat_non, group.by.vars = "orig.ident", dims.use = 1:17)

# Find clusters and run UMAP
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4)
##excitatory
seurat_ex <- FindNeighbors(object = seurat_ex, reduction = "harmony", dims=1:25)
seurat_ex <- FindClusters(object = seurat_ex, method="igraph", algorithm = 4, #use leiden algorithm
                                  resolution = resolutions)
seurat_ex <- RunUMAP(seurat_ex, reduction = "harmony", dims = 1:25)
##inhibitory
seurat_in <- FindNeighbors(object = seurat_in, reduction = "harmony", dims=1:20)
seurat_in <- FindClusters(object = seurat_in, method="igraph", algorithm = 4, #use leiden algorithm
                          resolution = resolutions)
seurat_in <- RunUMAP(seurat_in, reduction = "harmony", dims = 1:20)
##non-neuronal
resolutions <- c(0.4, 0.6, 0.8, 1.0)
seurat_non <- FindNeighbors(object = seurat_non, reduction = "harmony", dims=1:17)
seurat_non <- FindClusters(object = seurat_non, method="igraph", algorithm = 4, #use leiden algorithm
                          resolution = resolutions)
seurat_non <- RunUMAP(seurat_non, reduction = "harmony", dims = 1:17)

## 3.1 Excitatory ####
# Plot UMAP for each resolution
umap_plots <- list()
# Loop through each resolution and store the plot
for (res in resolutions) {
  # Generate UMAP plot and store in the list
  umap_plot <- DimPlot(seurat_ex, reduction = "umap", 
                       group.by = paste0("SCT_snn_res.",res),
                       label = TRUE, label.size = 5, repel = F) + 
    ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
  umap_plots[[as.character(res)]] <- umap_plot
}
# Combine UMAP plots into one plot for comparison with adjusted layout
combined_umap_plot <- wrap_plots(umap_plots, ncol = 2) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust margins to increase space between plots
  theme(legend.position = "bottom")  # Move legend to the bottom
# Display the combined UMAP plot
print(combined_umap_plot)
# Check clustree to determine resolution
p <- clustree(seurat_ex, prefix = "SCT_snn_res.")
p
# Check clustering by sample
res = 1.2
DimPlot(seurat_ex, reduction = "umap", 
        group.by = paste0("SCT_snn_res.",res),
        split.by = "orig.ident", ncol = 4,
        label = TRUE, label.size = 5, repel = F) + 
  ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
# Use resolution = 1.2 for excitatory
Idents(object = seurat_ex) <- "SCT_snn_res.1.2"
# Known markers
markers <- list(
  neuronal = c("RBFOX3", "FOXP2", "SYNPR"),
  non_neuronal = c("GFAP", "OLIG2", "CSF1R", "CLDN5", "PTPRC"),
  excitatory_general = c("SLC17A6", "SLC17A7"),
  excitatory_upper = c("CUX2", "GLRA3", "LAMP5"),
  excitatory_middle = c("BHLHE22", "RBP4", "RORB"),
  excitatory_lower = c("NR4A2", "TLE4"),
  inhibitory_general = c("CCK", "GAD1", "GAD2", "SLC32A1"),
  inhibitory_subtypes = c("GAD1", "LAMP5", "VIP", "SST", "PVALB"),
  additional = c("TBR1", "NR4A2", "RBP4", "PAX6", "SOX6", "KCNK2"),
  Ch_IN = c("GAD1","UNC5B","RORA","TRPS1","NFIB"),
  L4_5_IT = c("RORB","GPR26","TMSB10"),
  L5_corticobulbar_tract = c("GRIK1","SLC24A3"),
  L5_Betz_VEN = c("POU3F1"),
  L6b_CT = c("SULF1","SRGAP1","SEMA5A","SYT6","ZFHX3"),
  L6_IT_Car3 = c("NR4A2","SYT6","NTNG2","POSTN"),
  L6_IT = c("ITGB8","SNED1","OPRK1")
)
# Annotate clusters with known markers and visualize
for (marker_group in names(markers)) {
  marker_plots <- lapply(markers[[marker_group]], function(marker) {
    p <- FeaturePlot(seurat_ex,
                     reduction = "umap",
                     features = marker,
                     order = TRUE,
                     min.cutoff = 'q10',
                     label = TRUE) +
      ggtitle(marker)
    return(p)
  })
  combined_plot <- wrap_plots(marker_plots, ncol = 3) +
    plot_annotation(title = paste(marker_group, "markers")) & 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  print(combined_plot)
}
# Assign cell types ####
res <- "1.2"
Idents(object = seurat_ex) <- paste0("SCT_snn_res.",res)
nclusters <- nlevels(seurat_ex@active.ident)
##general subtypes
seurat_ex <- RenameIdents(object = seurat_ex, 
                                  "2" = "L2-3", "4" = "L2-3", "8" = "L2-3", "9" = "L2-3", "11" = "L2-3", "12" = "L2-3", "6" = "L2-3", "7" = "L2-3", 
                                  "3" = "L4-5", "21" = "L4-5", "14" = "L4-5", "10" = "L4-5", "18" = "L4-5", "5" = "L4-5", "13" = "L4-5", "16" = "L4-5", "17" = "L4-5", "22" = "L4-5",
                                  "1" = "L6", "19" = "L6", "20" = "L6", "15" = "L6") 
seurat_ex$subtype_general <- seurat_ex@active.ident
DimPlot(object = seurat_ex, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = FALSE)
##detailed subtypes
Idents(object = seurat_ex) <- paste0("SCT_snn_res.",res)
seurat_ex <- RenameIdents(object = seurat_ex, 
                                  "2" = "L2-3 IT", "4" = "L2-3 IT", "8" = "L2-3 IT", "9" = "L2-3 IT", "11" = "L2-3 IT", "12" = "L2-3 IT", "6" = "L2-3 IT", "7" = "L2-3 IT",
                                  "3" = "L4 IT", "21" = "L4 IT", "14" = "L4 IT", 
                                  "10" = "L5 IT", "18" = "L5 IT", "5" = "L5 IT", "13" = "L5 IT", "16" = "L5 IT",
                                  "17" = "L5 corticobulbar tract",
                                  "22" = "L5 Bentz VEN",
                                  "20" = "L6b",
                                  "15" = "L6 CT", 
                                  "19" = "L6 IT Car3", 
                                  "1" = "L6 IT") 
seurat_ex$subtype_detailed <- seurat_ex@active.ident
DimPlot(object = seurat_ex, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = FALSE)
# Find cluster markers
Idents(object = seurat_ex) <- "SCT_snn_res.1.2"
markers <- FindAllMarkers(seurat_ex, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p <- DoHeatmap(seurat_ex, features = top10$gene) + NoLegend() + ggtitle("Res = 1.2")
p

## 3.2 Inhibitory ####
# Plot UMAP for each resolution
umap_plots <- list()
# Loop through each resolution and store the plot
for (res in resolutions) {
  # Generate UMAP plot and store in the list
  umap_plot <- DimPlot(seurat_in, reduction = "umap", 
                       group.by = paste0("SCT_snn_res.",res),
                       label = TRUE, label.size = 5, repel = F) + 
    ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
  umap_plots[[as.character(res)]] <- umap_plot
}
# Combine UMAP plots into one plot for comparison with adjusted layout
combined_umap_plot <- wrap_plots(umap_plots, ncol = 2) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust margins to increase space between plots
  theme(legend.position = "bottom")  # Move legend to the bottom
# Display the combined UMAP plot
print(combined_umap_plot)
# Check clustree to determine resolution
p <- clustree(seurat_in, prefix = "SCT_snn_res.")
p
# Check clustering by sample
res = 1.4
DimPlot(seurat_in, reduction = "umap", 
        group.by = paste0("SCT_snn_res.",res),
        split.by = "orig.ident", ncol = 4,
        label = TRUE, label.size = 5, repel = F) + 
  ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
# Use resolution = 1.4 for inhibitory
Idents(object = seurat_in) <- "SCT_snn_res.1.4"
# Known markers
markers <- list(
  neuronal = c("RBFOX3", "FOXP2", "SYNPR"),
  non_neuronal = c("GFAP", "OLIG2", "CSF1R", "CLDN5", "PTPRC"),
  excitatory_general = c("SLC17A6", "SLC17A7"),
  excitatory_upper = c("CUX2", "GLRA3", "LAMP5"),
  excitatory_middle = c("BHLHE22", "RBP4", "RORB"),
  excitatory_lower = c("NR4A2", "TLE4"),
  inhibitory_general = c("CCK", "GAD1", "GAD2", "SLC32A1"),
  inhibitory_subtypes = c("GAD1", "LAMP5", "VIP", "SST", "PVALB"),
  additional = c("TBR1", "NR4A2", "RBP4", "PAX6", "SOX6", "KCNK2"),
  Ch_IN = c("GAD1","UNC5B","RORA","TRPS1","NFIB"),
  L4_5_IT = c("RORB","GPR26","TMSB10"),
  L5_corticobulbar_tract = c("GRIK1","SLC24A3"),
  L5_Betz_VEN = c("POU3F1"),
  L6b_CT = c("SULF1","SRGAP1","SEMA5A","SYT6","ZFHX3"),
  L6_IT_Car3 = c("NR4A2","SYT6","NTNG2","POSTN"),
  L6_IT = c("ITGB8","SNED1","OPRK1") 
)
# Annotate clusters with known markers and visualize
for (marker_group in names(markers)) {
  marker_plots <- lapply(markers[[marker_group]], function(marker) {
    p <- FeaturePlot(seurat_in,
                     reduction = "umap",
                     features = marker,
                     order = TRUE,
                     min.cutoff = 'q10',
                     label = TRUE) +
      ggtitle(marker)
    return(p)
  })
  combined_plot <- wrap_plots(marker_plots, ncol = 3) +
    plot_annotation(title = paste(marker_group, "markers")) & 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  print(combined_plot)
}
# Assign cell types ####
res <- "1.4"
Idents(object = seurat_in) <- paste0("SCT_snn_res.",res)
nclusters <- nlevels(seurat_in@active.ident)
##general subtypes
seurat_in$subtype_general <- "IN"
##detailed subtypes
Idents(object = seurat_in) <- paste0("SCT_snn_res.",res)
seurat_in <- RenameIdents(object = seurat_in, 
                          "18" = "Ch IN", 
                          "1" = "PV+ IN", "10" = "PV+ IN", "6" = "PV+ IN", "19" = "PV+ IN", "7" = "PV+ IN",
                          "4" = "SST+ IN", "3" = "SST+ IN", "5" = "SST+ IN", "17" = "SST+ IN", "24" = "SST+ IN", "23" = "SST+ IN", 
                          "12" = "VIP+ IN", "16" = "VIP+ IN", "8" = "VIP+ IN", "13" = "VIP+ IN", "9" = "VIP+ IN", 
                          "11" = "LAMP5+ IN", "2" = "LAMP5+ IN", "22" = "LAMP5+ IN", 
                          "15" = "SST+ LAMP5+ IN",
                          "14" = "VIP+like IN", "20" = "VIP+like IN", "21" = "VIP+like IN") 
seurat_in$subtype_detailed <- seurat_in@active.ident
DimPlot(object = seurat_in, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = FALSE)
# Find cluster markers
Idents(object = seurat_in) <- "SCT_snn_res.1.4"
markers <- FindAllMarkers(seurat_in, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p <- DoHeatmap(seurat_in, features = top10$gene) + NoLegend() + ggtitle("Res = 1.4")
p

## 3.3 Non-neuronal ####
# Plot UMAP for each resolution
umap_plots <- list()
# Loop through each resolution and store the plot
for (res in resolutions) {
  # Generate UMAP plot and store in the list
  umap_plot <- DimPlot(seurat_non, reduction = "umap", 
                       group.by = paste0("SCT_snn_res.",res),
                       label = TRUE, label.size = 5, repel = F) + 
    ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
  umap_plots[[as.character(res)]] <- umap_plot
}
# Combine UMAP plots into one plot for comparison with adjusted layout
combined_umap_plot <- wrap_plots(umap_plots, ncol = 2) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust margins to increase space between plots
  theme(legend.position = "bottom")  # Move legend to the bottom
# Display the combined UMAP plot
print(combined_umap_plot)
# Check clustree to determine resolution
p <- clustree(seurat_non, prefix = "SCT_snn_res.")
p
# Check clustering by sample
res = 0.4
DimPlot(seurat_non, reduction = "umap", 
        group.by = paste0("SCT_snn_res.",res),
        split.by = "orig.ident", ncol = 4,
        label = TRUE, label.size = 5, repel = F) + 
  ggtitle(paste("Resolution:", res)) + theme(legend.position = "none")
# Use resolution = 0.4 for nonneuronal
Idents(object = seurat_non) <- "SCT_snn_res.0.4"
# Known markers
markers <- list(
  neuronal = c("RBFOX3", "FOXP2", "SYNPR"),
  non_neuronal = c("GFAP", "OLIG2", "CSF1R", "CLDN5", "PTPRC"),
  oligo = c("OLIG2", "PLP1", "PDGFRA")
)
# Annotate clusters with known markers and visualize
for (marker_group in names(markers)) {
  marker_plots <- lapply(markers[[marker_group]], function(marker) {
    p <- FeaturePlot(seurat_non,
                     reduction = "umap",
                     features = marker,
                     order = TRUE,
                     min.cutoff = 'q10',
                     label = TRUE) +
      ggtitle(marker)
    return(p)
  })
  combined_plot <- wrap_plots(marker_plots, ncol = 3) +
    plot_annotation(title = paste(marker_group, "markers")) & 
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  print(combined_plot)
}
# Assign cell types ####
res <- "0.4"
Idents(object = seurat_non) <- paste0("SCT_snn_res.",res)
nclusters <- nlevels(seurat_non@active.ident)
##general subtypes
seurat_non$subtype_general <- "Non-neuronal"
##detailed subtypes
Idents(object = seurat_non) <- paste0("SCT_snn_res.",res)
seurat_non <- RenameIdents(object = seurat_non, 
                          "2" = "Astrocyte", 
                          "3" = "Oligo",
                          "6" = "OPC",
                          "4" = "Endothelial", "7" = "Endothelial", 
                          "5" = "Microglia/Blood") 
seurat_non$subtype_detailed <- seurat_non@active.ident
DimPlot(object = seurat_non, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = FALSE)
# Find cluster markers
Idents(object = seurat_non) <- "SCT_snn_res.0.4"
markers <- FindAllMarkers(seurat_non, assay = "SCT", only.pos = TRUE, recorrect_umi = FALSE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p <- DoHeatmap(seurat_non, features = top10$gene) + NoLegend() + ggtitle("Res = 0.4")
p


# 4. Assign sub cell type back to integrated seurat ####
res <- "0.8"
Idents(object = seurat_integrated) <- paste0("SCT_snn_res.",res)
nclusters <- nlevels(seurat_integrated@active.ident)
df <- data.frame(rowname = c(rownames(seurat_ex@meta.data), rownames(seurat_in@meta.data), rownames(seurat_non@meta.data)),
                 subtype_general = c(as.character(seurat_ex$subtype_general), as.character(seurat_in$subtype_general), as.character(seurat_non$subtype_general)),
                 subtype_detailed = c(as.character(seurat_ex$subtype_detailed), as.character(seurat_in$subtype_detailed), as.character(seurat_non$subtype_detailed))
)
meta <- seurat_integrated@meta.data
meta$celltype <- as.character(meta$celltype)
meta$celltype[match(df$rowname, rownames(meta))] <- df$subtype_general
meta$celltype_detailed <- meta$celltype
meta$celltype_detailed[match(df$rowname, rownames(meta))] <- as.character(df$subtype_detailed)
meta$celltype <- as.factor(meta$celltype)
meta$celltype_detailed <- as.factor(meta$celltype_detailed)
meta$clinical[which(meta$clinical=="Normal(Ctrl)")] <- "Control"
meta$brain_region[which(meta$brain_region=="BA9(PFC)")] <- "PFC"
meta$clinical_TDP43 <- paste(meta$clinical, meta$TDP43_status, sep="_")
meta$region_clinical_TDP43 <- paste(meta$brain_region, meta$clinical, meta$TDP43_status, sep="_")
seurat_integrated@meta.data <- meta
head(seurat_integrated@meta.data)

seurat_integrated@active.ident <- seurat_integrated$celltype
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = F)

seurat_integrated@active.ident <- seurat_integrated$celltype_detailed
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = F)


# 5. Extract neuronal cell fraction ####
df.plot <- FetchData(seurat_integrated, 
                     vars = c("orig.ident","case_ID","clinical","brain_region","TDP43_status","clinical_TDP43","region_clinical_TDP43","celltype_general","celltype","celltype_detailed"))
# Remove non-neuronal cells
df.plot$celltype_general <- as.character(df.plot$celltype_general)
df.plot <- df.plot[df.plot$celltype_general %in% c("Excitatory","Inhibitory"),]
# Count number of cells in each cell type
n_cells <- df.plot %>% group_by(orig.ident,case_ID,clinical,brain_region,TDP43_status,clinical_TDP43,region_clinical_TDP43,celltype) %>%
  count(celltype_detailed) %>% group_by(orig.ident) %>% mutate(n_by_sample=sum(n))
#celltype
n_cells <- n_cells %>% group_by(orig.ident,case_ID,clinical,brain_region,TDP43_status,clinical_TDP43,region_clinical_TDP43,celltype) %>%
  mutate(n_simple=sum(n)) %>% 
  group_by(celltype,region_clinical_TDP43) %>% mutate(frac_region_clinical_TDP43_simple=sum(n_simple)/sum(n_by_sample))
write.table(n_cells, "n_cells.neuronal.txt", sep="\t", quote=F, row.names=F)


# 6. Plot DEG of TOP1 genes ####
comparisons <- c("als_normal_vs_control","als_TDP_vs_control","ftd_normal_vs_control","ftd_TDP_vs_control","als_TDP_vs_als_normal","ftd_TDP_vs_ftd_normal")
target_genes <- c("RRM1", "RRM2B", "RNASEH2A", "RNASEH2B", "RNASEH2C")
nebula <- readRDS("deg_output.FTD_vs_control.nebula.rds")
celltypes <- names(nebula)
df2 <- c()
for (celltype in celltypes){
  tmp <- nebula[[celltype]]$summary
  tmp$p_val_adj <- p.adjust(tmp$`p_case_statusNormal(Ctrl)`, method = "BH")
  tmp <- tmp[tmp$gene %in% target_genes,]
  tmp$avg_log2FC <- -tmp$`logFC_case_statusNormal(Ctrl)` #convert to disease vs. control
  tmp <- tmp[,c("gene","avg_log2FC","p_val_adj")]
  tmp$celltype <- celltype
  df2 <- rbind(df2, tmp)
}
df2$comparison <- "FTD_vs_Control"
df <- rbind(df1, df2)
# Plot heatmap
df$comparison <- factor(df$comparison, levels = c("ALS_vs_Control","FTD_vs_Control"))
celltypes <- sort(unique(df$celltype))
celltypes <- celltypes[c(4:length(celltypes),2)] #neuronal celltypes only
plist <- list()
for (i in 1:length(unique(df$comparison))){
  mycomparison <- unique(df$comparison)[i]
  df.plot <- df[df$comparison==mycomparison&df$gene%in%target_genes&df$celltype%in%celltypes,c("gene","avg_log2FC","p_val_adj","celltype")]
  df.plot$mylabel <- ifelse(df.plot$p_val_adj<0.05, "*", "")
  df.plot$gene <- factor(df.plot$gene, levels = rev(target_genes))
  df.plot$celltype <- factor(df.plot$celltype, levels = celltypes)
  p <- ggplot(df.plot, aes(x=celltype, y=gene, fill=avg_log2FC)) +
    geom_tile() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    geom_text(aes(label=mylabel)) +
    labs(title=paste0(mycomparison, " (Nebula)")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
  plist[[i]] <- p
}
p <- grid.arrange(grobs = plist, ncol = 2)
p



