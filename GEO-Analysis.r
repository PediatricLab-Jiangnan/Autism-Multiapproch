# Load required packages
library(GEOquery)
library(limma)
library(sva)             # For batch correction
library(dplyr)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)    # Human gene annotation database
library(annotate)        # For gene annotation
library(preprocessCore)  # For quantile normalization
library(tidyr)
library(tibble)
library(gridExtra)

# 1. Download dataset and separate platforms -----------------------------
gset <- getGEO("GSE18123", destdir = ".", getGPL = TRUE)

# Automatically detect multi-platform data
if (length(gset) > 1) {
  gpl570 <- gset[[which(sapply(gset, annotation) == "GPL570")]]
  gpl6244 <- gset[[which(sapply(gset, annotation) == "GPL6244")]]
} else {
  stop("No multi-platform data detected. Please check GEO record.")
}

# 2. Expression matrix preprocessing -------------------------------------
expr570 <- exprs(gpl570)
expr6244 <- exprs(gpl6244)

cat("=== Raw Data Summary ===\n")
cat("GPL570 matrix dim:", dim(expr570), "\n")
cat("GPL570 value range:", range(expr570), "\n")
cat("GPL6244 matrix dim:", dim(expr6244), "\n")
cat("GPL6244 value range:", range(expr6244), "\n")

# 3. Get platform annotation info
gpl570_anno <- fData(gpl570)
gpl6244_anno <- fData(gpl6244)

# 4. Gene annotation functions -------------------------------------------
# Improved GPL570 annotation function
annotate_gpl570 <- function(expr_matrix, platform_anno) {
  cat("Annotating GPL570 using Gene Symbol column\n")
  gene_symbols <- platform_anno$`Gene Symbol`
  names(gene_symbols) <- rownames(platform_anno)
  valid_probes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---" & !grepl("///", gene_symbols)
  expr_filtered <- expr_matrix[valid_probes, ]
  gene_symbols_filtered <- gene_symbols[valid_probes]
  expr_df <- data.frame(expr_filtered, check.names = FALSE)
  expr_df$gene_symbol <- gene_symbols_filtered
  expr_collapsed <- expr_df %>%
    group_by(gene_symbol) %>%
    summarise_all(median, na.rm = TRUE) %>%
    column_to_rownames("gene_symbol")
  return(as.matrix(expr_collapsed))
}

# Improved GPL6244 annotation function
annotate_gpl6244 <- function(expr_matrix, platform_anno) {
  cat("Annotating GPL6244 using gene_assignment column\n")
  gene_assignments <- platform_anno$gene_assignment
  extract_gene_symbol <- function(assignment_string) {
    if(is.na(assignment_string) || assignment_string == "" || assignment_string == "---") return(NA)
    genes <- unlist(strsplit(assignment_string, " /// "))
    first_gene <- genes[1]
    gene_parts <- unlist(strsplit(first_gene, " // "))
    if(length(gene_parts) >= 2) return(trimws(gene_parts[2])) else return(NA)
  }
  gene_symbols <- sapply(gene_assignments, extract_gene_symbol)
  names(gene_symbols) <- rownames(platform_anno)
  valid_probes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  # Try alternative annotation if too few valid probes
  if(sum(valid_probes) < 1000) {
    probe_ids <- rownames(platform_anno)
    tryCatch({
      library(illuminaHumanv4.db)
      gene_symbols_alt <- getSYMBOL(probe_ids, data = 'illuminaHumanv4')
      valid_alt <- !is.na(gene_symbols_alt)
      gene_symbols[valid_alt] <- gene_symbols_alt[valid_alt]
    }, error = function(e) {
      tryCatch({
        gene_symbols_alt <- getSYMBOL(probe_ids, data = 'org.Hs.eg')
        valid_alt <- !is.na(gene_symbols_alt)
        gene_symbols[valid_alt] <- gene_symbols_alt[valid_alt]
      }, error = function(e2) {})
    })
    valid_probes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  }
  expr_filtered <- expr_matrix[valid_probes, ]
  gene_symbols_filtered <- gene_symbols[valid_probes]
  expr_df <- data.frame(expr_filtered, check.names = FALSE)
  expr_df$gene_symbol <- gene_symbols_filtered
  expr_collapsed <- expr_df %>%
    group_by(gene_symbol) %>%
    summarise_all(median, na.rm = TRUE) %>%
    column_to_rownames("gene_symbol")
  return(as.matrix(expr_collapsed))
}

# Annotate gene symbols for both platforms
cat("\nStarting GPL570 annotation...\n")
expr570_annotated <- annotate_gpl570(expr570, gpl570_anno)
cat("\nStarting GPL6244 annotation...\n")
expr6244_annotated <- annotate_gpl6244(expr6244, gpl6244_anno)

# 5. Expression normalization --------------------------------------------
cat("\n=== Normalization ===\n")
normalize_expression <- function(expr_matrix, method = "quantile", platform_name) {
  cat("Normalizing", platform_name, "using", method, "...\n")
  if(method == "quantile") {
    expr_normalized <- normalize.quantiles(expr_matrix)
    rownames(expr_normalized) <- rownames(expr_matrix)
    colnames(expr_normalized) <- colnames(expr_matrix)
  } else if(method == "log2") {
    data_range <- max(expr_matrix, na.rm = TRUE) - min(expr_matrix, na.rm = TRUE)
    if(data_range > 50) {
      expr_normalized <- log2(expr_matrix + 1)
    } else {
      expr_normalized <- expr_matrix
    }
  } else if(method == "zscore") {
    expr_normalized <- t(scale(t(expr_matrix)))
  } else if(method == "robust") {
    expr_normalized <- expr_matrix
    for(i in 1:nrow(expr_matrix)) {
      gene_values <- expr_matrix[i, ]
      median_val <- median(gene_values, na.rm = TRUE)
      mad_val <- mad(gene_values, na.rm = TRUE)
      if(mad_val > 0) expr_normalized[i, ] <- (gene_values - median_val) / mad_val
    }
  }
  return(expr_normalized)
}

normalization_method <- "log2"  # Choose from: quantile, log2, zscore, robust
expr570_normalized <- normalize_expression(expr570_annotated, method = normalization_method, platform_name = "GPL570")
expr6244_normalized <- normalize_expression(expr6244_annotated, method = normalization_method, platform_name = "GPL6244")

# 6. Visualization of normalization effect -------------------------------
create_normalization_comparison <- function(expr_before, expr_after, platform_name) {
  set.seed(123)
  n_samples <- min(20, ncol(expr_before))
  selected_samples <- sample(colnames(expr_before), n_samples)
  before_data <- expr_before[, selected_samples] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(status = "Before")
  after_data <- expr_after[, selected_samples] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(status = "After")
  plot_data <- rbind(before_data, after_data)
  plot_data$status <- factor(plot_data$status, levels = c("Before", "After"))
  p <- ggplot(plot_data, aes(x = sample, y = expression, fill = status)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste0(platform_name, " - Normalization Comparison"),
         x = "Sample", y = "Expression", fill = "Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    facet_wrap(~status, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Before" = "#FF6B6B", "After" = "#4ECDC4"))
  return(p)
}
norm_comparison_570 <- create_normalization_comparison(expr570_annotated, expr570_normalized, "GPL570")
norm_comparison_6244 <- create_normalization_comparison(expr6244_annotated, expr6244_normalized, "GPL6244")
print(norm_comparison_570)
print(norm_comparison_6244)
ggsave("GPL570_normalization_comparison.pdf", norm_comparison_570, width = 12, height = 8, dpi = 300)
ggsave("GPL6244_normalization_comparison.pdf", norm_comparison_6244, width = 12, height = 8, dpi = 300)

# 7. Get common genes ----------------------------------------------------
common_genes <- intersect(rownames(expr570_normalized), rownames(expr6244_normalized))
cat("\nNumber of common genes:", length(common_genes), "\n")
if(length(common_genes) < 5000) {
  cat("\nWarning: Few common genes (<5000). Check annotation and platform compatibility.\n")
}
expr570_common <- expr570_normalized[common_genes, ]
expr6244_common <- expr6244_normalized[common_genes, ]

# 8. Extract phenotype and group -----------------------------------------
gpl570_pdata <- pData(gpl570)
gpl570_pdata$group <- ifelse(grepl("diagnosis: CONTROL", gpl570_pdata$characteristics_ch1), "CONTROL", "AUTISM")
gpl570_pdata$group <- factor(gpl570_pdata$group, levels = c("CONTROL", "AUTISM"))
gpl6244_pdata <- pData(gpl6244)
gpl6244_pdata$group <- ifelse(grepl("diagnosis: CONTROL", gpl6244_pdata$characteristics_ch1), "CONTROL", "AUTISM")
gpl6244_pdata$group <- factor(gpl6244_pdata$group, levels = c("CONTROL", "AUTISM"))

# 9. Merge datasets ------------------------------------------------------
expr570_common <- expr570_common[, rownames(gpl570_pdata)]
expr6244_common <- expr6244_common[, rownames(gpl6244_pdata)]
combined_expr <- cbind(expr570_common, expr6244_common)
gpl570_pdata$platform <- "GPL570"; gpl6244_pdata$platform <- "GPL6244"
gpl570_pdata$batch <- 1; gpl6244_pdata$batch <- 2
common_cols <- intersect(colnames(gpl570_pdata), colnames(gpl6244_pdata))
gpl570_pdata_subset <- gpl570_pdata[, common_cols]
gpl6244_pdata_subset <- gpl6244_pdata[, common_cols]
combined_pdata <- rbind(gpl570_pdata_subset, gpl6244_pdata_subset)
if(!all(colnames(combined_expr) == rownames(combined_pdata))) {
  combined_pdata <- combined_pdata[colnames(combined_expr), ]
}

# 10. PCA visualization before batch correction --------------------------
perform_pca <- function(expr_matrix, pdata, title_prefix) {
  pca_data <- t(expr_matrix)
  pca_result <- prcomp(pca_data, scale. = TRUE)
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  pca_df <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    group = pdata$group,
    platform = pdata$platform,
    sample_id = rownames(pdata)
  )
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = platform)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = paste0(title_prefix, " - PCA"),
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = platform, shape = group)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = paste0(title_prefix, " - Color by Platform"),
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  return(list(p1 = p1, p2 = p2, pca_result = pca_result))
}
pca_before <- perform_pca(combined_expr, combined_pdata, "Before Batch Correction")
print(pca_before$p1)
print(pca_before$p2)

# 11. Batch correction ---------------------------------------------------
batch <- combined_pdata$batch
group <- combined_pdata$group
mod <- model.matrix(~group, data = combined_pdata)
combined_expr_corrected <- ComBat(dat = combined_expr, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
cat("Batch correction completed\n")

# 12. PCA visualization after batch correction ---------------------------
pca_after <- perform_pca(combined_expr_corrected, combined_pdata, "After Batch Correction")
print(pca_after$p1)
print(pca_after$p2)
pca_comparison <- grid.arrange(
  pca_before$p2 + ggtitle("Before Batch Correction - Color by Platform"),
  pca_after$p2 + ggtitle("After Batch Correction - Color by Platform"),
  ncol = 2
)
ggsave("PCA_batch_correction_comparison.png", pca_comparison, width = 12, height = 6, dpi = 300)

# 13. Export data --------------------------------------------------------
write.csv(combined_expr, "GSE18123_expression_before_correction.csv", row.names = TRUE)
write.csv(combined_expr_corrected, "GSE18123_expression_after_correction.csv", row.names = TRUE)
write.csv(combined_pdata, "GSE18123_sample_info.csv", row.names = TRUE)

# 14. Detailed PCA analysis with ellipse ---------------------------------
perform_detailed_pca <- function(expr_matrix, pdata) {
  pca_data <- t(expr_matrix)
  pca_result <- prcomp(pca_data, scale. = TRUE)
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  pca_df <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    PC3 = pca_result$x[,3],
    group = pdata$group,
    platform = pdata$platform,
    sample_id = rownames(pdata)
  )
  return(list(pca_df = pca_df, var_explained = var_explained, pca_result = pca_result))
}
detailed_pca <- perform_detailed_pca(combined_expr_corrected, combined_pdata)
pca_df <- detailed_pca$pca_df
var_explained <- detailed_pca$var_explained
pca_group_ellipse <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = group), alpha = 0.2, level = 0.95, geom = "polygon") +
  scale_color_manual(values = c("CONTROL" = "#3498DB", "AUTISM" = "#E74C3C")) +
  scale_fill_manual(values = c("CONTROL" = "#3498DB", "AUTISM" = "#E74C3C")) +
  labs(
    title = "PCA - After Batch Correction (By Diagnosis)",
    subtitle = "95% Confidence Ellipse",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(title = "Diagnosis", override.aes = list(size = 4)),
    fill = guide_legend(title = "Diagnosis")
  )
print(pca_group_ellipse)
