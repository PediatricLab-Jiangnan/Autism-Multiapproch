# Load required packages
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(hdf5r)   

# Function for initial QC and normalization for a single sample
process_single_sample <- function(h5_file_path, sample_name, condition) {
  # Read data from h5 file and create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = Read10X_h5(h5_file_path),
    project = sample_name
  )
  # Add metadata
  seurat_obj$batch <- sample_name
  seurat_obj$condition <- condition
  # Calculate mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
}

# Read and process each sample
sample_list <- list()
for(i in 1:3) {
  # Control samples
  con_name <- paste0("Con", i)
  con_file <- paste0(con_name, ".h5")
  sample_list[[con_name]] <- process_single_sample(con_file, con_name, "Control")
  # Autism samples
  se_name <- paste0("ASD", i)
  se_file <- paste0(se_name, ".h5")
  sample_list[[se_name]] <- process_single_sample(se_file, se_name, "Autism")
}

# Merge all samples for visualization
merged_seurat <- merge(sample_list[[1]], 
                       y = sample_list[2:length(sample_list)], 
                       add.cell.ids = names(sample_list))

# Pre-QC violin plots
qc_plots <- VlnPlot(merged_seurat, 
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "batch",
                    ncol = 3,
                    pt.size = 0) 
qc_plots[[1]] <- qc_plots[[1]] + ggtitle("nFeature_RNA")
qc_plots[[2]] <- qc_plots[[2]] + ggtitle("nCount_RNA")
qc_plots[[3]] <- qc_plots[[3]] + ggtitle("percent.mt")
combined_plot <- wrap_plots(qc_plots)
pdf("pre_qc_violins_combined.pdf", width = 15, height = 6)
print(combined_plot)
dev.off()

# Step 2: Quality Control
for(sample_name in names(sample_list)) {
  sample_list[[sample_name]] <- subset(sample_list[[sample_name]], 
                                       subset = nFeature_RNA > 200 & 
                                         nFeature_RNA < 6000 & 
                                         percent.mt < 20)
}

# Merge QCed data
merged_seurat <- merge(sample_list[[1]], 
                       y = sample_list[2:length(sample_list)], 
                       add.cell.ids = names(sample_list))

# Post-QC violin plots
post_qc_plots <- VlnPlot(merged_seurat, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         group.by = "batch",
                         ncol = 3,
                         pt.size = 0) 
post_qc_plots[[1]] <- post_qc_plots[[1]] + ggtitle("nFeature_RNA")
post_qc_plots[[2]] <- post_qc_plots[[2]] + ggtitle("nCount_RNA")
post_qc_plots[[3]] <- post_qc_plots[[3]] + ggtitle("percent.mt")
combined_post_qc_plot <- wrap_plots(post_qc_plots)
pdf("post_qc_violins_combined.pdf", width = 15, height = 6)
print(combined_post_qc_plot)
dev.off()

# Clean up
rm(sample_list, post_qc_plots, combined_post_qc_plot)
gc()

# Normalization and feature selection
merged_seurat <- merged_seurat |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData(vars.to.regress = c("percent.mt")) # Regress out mitochondrial effect

# PCA
merged_seurat <- RunPCA(merged_seurat, npcs = 30)

# Batch correction using Harmony
merged_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

# Save integrated Seurat object
saveRDS(merged_seurat, file = "merged_seurat_harmony.rds")

# Elbow plot for PCA
pdf("pca_elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(merged_seurat, reduction = "pca", ndims = 30) +
  ggtitle("PCA Elbow Plot") +
  theme_classic()
dev.off()

# Step 3: Clustering with different resolutions
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:20)
for(res in seq(0.1, 1.0, by = 0.1)) {
  merged_seurat <- FindClusters(merged_seurat, resolution = res, algorithm = 1)
}
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:20)

# Clustree analysis
library(clustree)
pdf("clustree_analysis.pdf", width = 12, height = 16)
clustree(merged_seurat, prefix = "RNA_snn_res.")
dev.off()

# UMAP plots for different resolutions
plot_list <- list()
for(res in seq(0.1, 1.0, by = 0.1)) {
  p <- DimPlot(merged_seurat, 
               reduction = "umap", 
               group.by = paste0("RNA_snn_res.", res),
               label = TRUE) +
    ggtitle(paste("Resolution:", res))
  plot_list[[length(plot_list) + 1]] <- p
}
pdf("umap_all_resolutions.pdf", width = 20, height = 15)
wrap_plots(plot_list, ncol = 3)
dev.off()

# Print cluster numbers for each resolution
cat("Number of clusters at each resolution:\n")
for(res in seq(0.1, 1.0, by = 0.1)) {
  clusters <- length(unique(merged_seurat[[paste0("RNA_snn_res.", res)]][,1]))
  cat(sprintf("Resolution %.1f: %d clusters\n", res, clusters))
}

# Cell annotation
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.1)
Idents(merged_seurat) <- "RNA_snn_res.0.1"

# Find marker genes for each cluster
markers <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Filter significant markers and sort
markers_sorted <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(markers, "all_markers_res0.1.csv", row.names = FALSE)
write.csv(markers_sorted, "top10_markers_res0.1.csv", row.names = FALSE)

# Heatmap for top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene)
pdf("markers_heatmap_res0.1.pdf", width = 12, height = 8)
DoHeatmap(merged_seurat, features = top_markers) +
  theme(axis.text.y = element_text(size = 8))
dev.off()

# Violin plot for top marker genes
pdf("markers_vlnplot_res0.1.pdf", width = 15, height = 10)
VlnPlot(merged_seurat, 
        features = unique(markers_sorted$gene[1:10]), 
        ncol = 5)
dev.off()

# Cluster size output
print(table(Idents(merged_seurat)))

# Automated cell type annotation based on marker genes
cell_type_markers <- list(
  "Oligodendrocytes" = c("OPALIN", "MAG", "MBP", "PLP1"),
  "Neurons" = c("FAM163A", "LINC00508", "TSHZ2", "IL1RAPL2", "EBF1", "LEF1", "TRABD2A", "HS3ST2", "SEMA3E", "SULF1", "SMYD1", "THEMIS"),
  "Astrocytes" = c("ETNPPL", "GJA1", "AQP4", "GFAP", "S100B"),
  "Microglia" = c("P2RY12", "CSF1R", "CX3CR1", "TMEM119"),
  "Excitatory_Neuron" = c("SLC17A7", "NEFM", "CAMK2A", "GRIN1"),
  "Oligodendrocyte_precursor_cells" = c("PDGFRA", "VCAN", "NG2", "OLIG2"),
  "Inhibitory_Neuron" = c("SST", "LHX6", "VIP", "CALB2", "KIT", "SV2C", "GAD1", "GAD2")
)

# Check marker genes present in data
all_markers <- unique(unlist(cell_type_markers))
present_markers <- all_markers[all_markers %in% rownames(merged_seurat)]
missing_markers <- all_markers[!all_markers %in% rownames(merged_seurat)]
cat("Total marker genes:", length(all_markers), "\n")
cat("Present marker genes:", length(present_markers), "\n")
cat("Missing marker genes:", length(missing_markers), "\n")

# Filter marker list to only present genes
cell_type_markers_filtered <- lapply(cell_type_markers, function(genes) {
  genes[genes %in% rownames(merged_seurat)]
})
cell_type_markers_filtered <- cell_type_markers_filtered[sapply(cell_type_markers_filtered, length) > 0]

# Average expression of marker genes per cluster
cluster_avg_exp <- AverageExpression(merged_seurat, 
                                     features = present_markers,
                                     group.by = "RNA_snn_res.0.1")$RNA

# Compute cell type scores for each cluster
cluster_scores <- data.frame(cluster = colnames(cluster_avg_exp))
rownames(cluster_scores) <- cluster_scores$cluster
for(cell_type in names(cell_type_markers_filtered)) {
  markers <- cell_type_markers_filtered[[cell_type]]
  if(length(markers) > 0) {
    scores <- if(length(markers) == 1) cluster_avg_exp[markers, ] else colMeans(cluster_avg_exp[markers, , drop = FALSE])
    cluster_scores[[cell_type]] <- scores
  }
}

# Assign cell type with highest score to each cluster
cluster_annotations <- apply(cluster_scores[, -1, drop = FALSE], 1, function(x) {
  if(all(is.na(x))) return("Unknown")
  cell_type <- names(which.max(x))
  return(cell_type)
})

# Map clusters to cell types, ensure names match
cluster_to_celltype <- cluster_annotations
names(cluster_to_celltype) <- gsub("^g", "", names(cluster_to_celltype))
current_clusters <- levels(Idents(merged_seurat))
mapping_clusters <- names(cluster_to_celltype)
missing_in_mapping <- setdiff(current_clusters, mapping_clusters)
extra_in_mapping <- setdiff(mapping_clusters, current_clusters)
if(length(missing_in_mapping) > 0) {
  for(missing_cluster in missing_in_mapping) {
    cluster_to_celltype[missing_cluster] <- "Unknown"
  }
}
if(length(extra_in_mapping) > 0) {
  cluster_to_celltype <- cluster_to_celltype[names(cluster_to_celltype) %in% current_clusters]
}

# Apply annotation
merged_seurat <- RenameIdents(merged_seurat, cluster_to_celltype)
merged_seurat$cell_type_auto <- Idents(merged_seurat)

# Summary statistics
cell_type_counts <- table(merged_seurat$cell_type_auto)
mapping_table <- table(merged_seurat$RNA_snn_res.0.1, merged_seurat$cell_type_auto)
annotation_summary <- data.frame(
  Cluster = names(cluster_to_celltype),
  Cell_Type = cluster_to_celltype,
  Cell_Count = as.numeric(table(merged_seurat$RNA_snn_res.0.1)[names(cluster_to_celltype)]),
  stringsAsFactors = FALSE
)
cell_type_summary <- aggregate(Cell_Count ~ Cell_Type, data = annotation_summary, sum)
cell_type_summary <- cell_type_summary[order(cell_type_summary$Cell_Count, decreasing = TRUE), ]

# Visualization of annotation
pdf("umap_auto_annotation.pdf", width = 16, height = 6)
p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE, label.size = 3, pt.size = 0.5) +
  ggtitle("Original Clusters (Resolution 0.1)") + theme_classic()
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type_auto", label = TRUE, label.size = 3, pt.size = 0.5) +
  ggtitle("Auto-annotated Cell Types") + theme_classic()
print(p1 + p2)
dev.off()

# Heatmap for cell type scores
library(pheatmap)
score_matrix <- as.matrix(cluster_scores[, -1])
rownames(score_matrix) <- gsub("^g", "", rownames(score_matrix))
pdf("cluster_celltype_scores_heatmap.pdf", width = 12, height = 8)
pheatmap(t(score_matrix), cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Cell Type Scores for Each Cluster", fontsize = 8, cellwidth = 15, cellheight = 15)
dev.off()

# Marker gene validation plots
pdf("marker_validation.pdf", width = 15, height = 20)
library(patchwork)
for(cell_type in names(cell_type_markers_filtered)) {
  markers <- cell_type_markers_filtered[[cell_type]]
  if(length(markers) > 0) {
    markers_to_plot <- markers[1:min(4, length(markers))]
    p <- FeaturePlot(merged_seurat, features = markers_to_plot, ncol = 2, pt.size = 0.5) +
      plot_annotation(title = paste("Markers for", gsub("_", " ", cell_type)))
    print(p)
  }
}
dev.off()

# Violin plots for marker genes by cell type
pdf("celltype_markers_violin.pdf", width = 15, height = 10)
for(cell_type in names(cell_type_markers_filtered)) {
  markers <- cell_type_markers_filtered[[cell_type]]
  if(length(markers) > 0) {
    markers_to_plot <- markers[1:min(2, length(markers))]
    p <- VlnPlot(merged_seurat, features = markers_to_plot, group.by = "cell_type_auto", ncol = 1, pt.size = 0) +
      plot_annotation(title = paste("Expression of", gsub("_", " ", cell_type), "markers"))
    print(p)
  }
}
dev.off()

# Save results
saveRDS(merged_seurat, "merged_seurat_auto_annotated.rds")
write.csv(cluster_scores, "cluster_celltype_scores.csv", row.names = TRUE)
write.csv(annotation_summary, "cluster_annotation_summary.csv", row.names = FALSE)
write.csv(cell_type_summary, "celltype_summary.csv", row.names = FALSE)
write.csv(mapping_table, "cluster_to_celltype_mapping.csv", row.names = TRUE)

# Final summary
cat("=== Auto annotation summary ===\n")
cat("Original clusters:", length(current_clusters), "\n")
cat("Annotated cell types:", length(unique(merged_seurat$cell_type_auto)), "\n")
cat("Cell type distribution:\n")
print(sort(table(merged_seurat$cell_type_auto), decreasing = TRUE))
cat("All result files saved.\n")

# Cell type distribution
print("Number of cells in each cell type:")
print(table(merged_seurat$cell_type_auto))

# DotPlot for marker genes
dotplot_pretty <- DotPlot(merged_seurat, features = present_markers, group.by = "cell_type_auto") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3, title = "Percent Expressed")) +
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colours = c('#330066', '#336699', '#66CC66', '#FFCC33'),
                        name = "Average Expression") +
  ggtitle("Cell Type Marker Genes Expression")
ggsave("dotplot_pretty.pdf", dotplot_pretty, width = 12, height = 10)

# UMAP visualizations
p1 <- DimPlot(merged_seurat, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  ggtitle("Clusters in UMAP") +
  xlab("UMAP_1") + ylab("UMAP_2")
ggsave("umap_clusters.pdf", p1, width = 10, height = 8)

p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "condition") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  ggtitle("Control vs Autism") +
  xlab("UMAP_1") + ylab("UMAP_2")
ggsave("umap_condition.pdf", p2, width = 10, height = 8)

p3 <- DimPlot(merged_seurat, reduction = "umap", group.by = "batch") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  ggtitle("Samples in UMAP") +
  xlab("UMAP_1") + ylab("UMAP_2")
ggsave("umap_batch.pdf", p3, width = 10, height = 8)

merged_seurat$condition <- factor(merged_seurat$condition, levels = c("Control", "Autism"))
p4 <- DimPlot(merged_seurat, reduction = "umap", split.by = "condition", label = TRUE, label.size = 4, repel = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  ggtitle("Clusters Split by Condition")
ggsave("umap_split_by_condition.pdf", p4, width = 15, height = 7)
