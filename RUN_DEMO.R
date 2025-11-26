# ============================================================================
# DEMONSTRATION: scRNA-seq Analysis Pipeline
# ============================================================================
# This script demonstrates the analysis workflow
# Real data files are in data/raw/ ready for analysis
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat("  📊 scRNA-seq ANALYSIS PIPELINE DEMONSTRATION\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat("\n")

# Set working directory
setwd("C:/Users/krupa/Desktop/Shiny_project")

# Load required libraries
library(tidyverse)
library(data.table)

# Define rowVars if not available
if (!exists("rowVars")) {
  rowVars <- function(x) {
    rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
  }
}

# ============================================================================
# STEP 1: Data Overview
# ============================================================================
cat("STEP 1: Available Data Files\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

data_files <- list.files("data/raw", pattern = "\\.gz$", full.names = FALSE)
datasets <- list(
  GSE123813 = list(
    name = "Basal Cell Carcinoma (BCC)",
    files = grep("GSE123813", data_files, value = TRUE),
    desc = "T-cell and immune microenvironment from BCC tumors"
  ),
  GSE243013 = list(
    name = "NSCLC Immune Microenvironment",
    files = grep("GSE243013", data_files, value = TRUE),
    desc = "Immune cells from non-small cell lung cancer"
  )
)

for (dataset in names(datasets)) {
  cat(sprintf("✓ %s: %s\n", dataset, datasets[[dataset]]$name))
  cat(sprintf("  Description: %s\n", datasets[[dataset]]$desc))
  cat(sprintf("  Files: %d\n", length(datasets[[dataset]]$files)))
  for (f in datasets[[dataset]]$files) {
    file_size <- file.size(paste0("data/raw/", f)) / (1024^2)
    cat(sprintf("    - %s (%.1f MB)\n", f, file_size))
  }
  cat("\n")
}

# ============================================================================
# STEP 2: Simulate scRNA-seq Data
# ============================================================================
cat("STEP 2: Creating Demonstration Dataset\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

set.seed(42)

# Parameters
n_cells <- 1000
n_genes <- 2000
n_cell_types <- 5

cat(sprintf("Simulating %d cells × %d genes...\n", n_cells, n_genes))

# Gene expression matrix
expr_matrix <- matrix(
  rpois(n_cells * n_genes, lambda = 5),
  nrow = n_genes,
  ncol = n_cells,
  dimnames = list(
    paste0("gene_", 1:n_genes),
    paste0("cell_", 1:n_cells)
  )
)

# Cell metadata
cell_types <- c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial")
metadata <- data.frame(
  cell_id = colnames(expr_matrix),
  n_counts = colSums(expr_matrix),
  n_genes = colSums(expr_matrix > 0),
  cell_type = sample(cell_types, n_cells, replace = TRUE),
  batch = sample(c("Batch1", "Batch2"), n_cells, replace = TRUE),
  stringsAsFactors = FALSE
)

cat(sprintf("✓ Expression matrix: %d genes × %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("✓ Sparsity: %.1f%%\n", 100 * sum(expr_matrix == 0) / length(expr_matrix)))

# ============================================================================
# STEP 3: Quality Control
# ============================================================================
cat("\nSTEP 3: Quality Control Metrics\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

qc_stats <- metadata %>%
  group_by(cell_type) %>%
  summarise(
    n_cells = n(),
    mean_counts = round(mean(n_counts), 1),
    median_genes = median(n_genes),
    max_genes = max(n_genes),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_cells))

cat("Cell Type Summary:\n")
print(as.data.frame(qc_stats), row.names = FALSE)

# Quality metrics summary
total_counts <- sum(expr_matrix)
total_genes <- sum(rowSums(expr_matrix) > 0)

cat(sprintf("\nGlobal Metrics:\n"))
cat(sprintf("  Total counts: %d\n", total_counts))
cat(sprintf("  Genes detected: %d\n", total_genes))
cat(sprintf("  Mean counts/cell: %.1f\n", mean(metadata$n_counts)))
cat(sprintf("  Mean genes/cell: %.1f\n", mean(metadata$n_genes)))

# ============================================================================
# STEP 4: Data Normalization & HVG Selection
# ============================================================================
cat("\nSTEP 4: Data Normalization\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

# Normalize (log-CPM)
cat("  Normalizing counts (CPM + log2)...\n")
expr_norm <- sweep(expr_matrix, 2, colSums(expr_matrix), "/") * 1e6
expr_log <- log2(expr_norm + 1)

# Calculate variance
gene_var <- rowVars(expr_log)
gene_mean <- rowMeans(expr_log)

# Select highly variable genes
n_hvg <- min(500, nrow(expr_log))
hvg_idx <- order(gene_var, decreasing = TRUE)[1:n_hvg]

cat(sprintf("✓ Selected %d highly variable genes\n", n_hvg))
cat(sprintf("  Mean expression range: [%.2f, %.2f]\n", 
            min(gene_mean[hvg_idx]), max(gene_mean[hvg_idx])))

# ============================================================================
# STEP 5: Dimensionality Reduction (PCA)
# ============================================================================
cat("\nSTEP 5: Dimensionality Reduction (PCA)\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

cat("  Computing PCA on top variable genes...\n")
pca_data <- t(expr_log[hvg_idx, ])
pca <- prcomp(pca_data, scale = TRUE, rank. = 30)

# Variance explained
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

cat(sprintf("✓ PCA computed (%d components)\n", length(pca$sdev)))
cat(sprintf("  PC1 explains: %.1f%%\n", var_explained[1]))
cat(sprintf("  PC2 explains: %.1f%%\n", var_explained[2]))
cat(sprintf("  Top 10 PCs explain: %.1f%%\n", sum(var_explained[1:10])))

# Create PCA dataframe for visualization
pca_df <- as.data.frame(pca$x[, 1:3])
pca_df$cell_type <- metadata$cell_type

# ============================================================================
# STEP 6: Clustering
# ============================================================================
cat("\nSTEP 6: Clustering Analysis\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

cat("  Running k-means clustering (k=5)...\n")
set.seed(42)
kmeans_result <- kmeans(pca$x[, 1:10], centers = 5, nstart = 10)

metadata$cluster <- as.factor(kmeans_result$cluster)
pca_df$cluster <- metadata$cluster

cluster_stats <- metadata %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    dominant_type = names(which.max(table(cell_type))),
    n_types = n_distinct(cell_type),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_cells))

cat("✓ Clustering complete:\n")
print(as.data.frame(cluster_stats), row.names = FALSE)

# ============================================================================
# STEP 7: Marker Gene Analysis
# ============================================================================
cat("\nSTEP 7: Marker Gene Identification\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

cat("  Computing marker genes per cluster...\n")

marker_list <- list()
for (clust in unique(metadata$cluster)) {
  in_cluster_idx <- which(metadata$cluster == clust)
  out_cluster_idx <- which(metadata$cluster != clust)
  
  if (length(in_cluster_idx) > 0 && length(out_cluster_idx) > 0) {
    # Mean expression
    mean_in <- rowMeans(expr_matrix[, in_cluster_idx, drop = FALSE])
    mean_out <- rowMeans(expr_matrix[, out_cluster_idx, drop = FALSE])
    
    # Log fold change
    lfc <- log2((mean_in + 1) / (mean_out + 1))
    
    # Get top markers
    top_idx <- order(lfc, decreasing = TRUE)[1:10]
    marker_list[[paste0("Cluster_", clust)]] <- rownames(expr_matrix)[top_idx]
  }
}

cat("✓ Top markers per cluster:\n")
for (clust_name in names(marker_list)) {
  markers <- marker_list[[clust_name]][1:5]
  cat(sprintf("  %s: %s\n", clust_name, paste(markers, collapse = ", ")))
}

# ============================================================================
# STEP 8: Save Results
# ============================================================================
cat("\nSTEP 8: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────────────\n\n")

# Save analysis object
analysis_result <- list(
  metadata = metadata,
  expression_matrix = expr_matrix,
  normalized_matrix = expr_log,
  hvg_indices = hvg_idx,
  pca = pca,
  pca_data = pca_df,
  clusters = metadata$cluster,
  cluster_stats = cluster_stats,
  marker_genes = marker_list
)

saveRDS(analysis_result, "data/processed/demo_analysis.rds")
cat("✓ Saved: data/processed/demo_analysis.rds\n")

# Save metadata
write.csv(metadata, "data/processed/cell_metadata.csv", row.names = FALSE)
cat("✓ Saved: data/processed/cell_metadata.csv\n")

# Save marker genes
marker_df <- data.frame(
  cluster = rep(names(marker_list), sapply(marker_list, length)),
  marker_gene = unlist(marker_list),
  stringsAsFactors = FALSE
)
write.csv(marker_df, "data/processed/marker_genes.csv", row.names = FALSE)
cat("✓ Saved: data/processed/marker_genes.csv\n")

# Save PCA coordinates
write.csv(pca_df, "data/processed/pca_coordinates.csv", row.names = FALSE)
cat("✓ Saved: data/processed/pca_coordinates.csv\n")

# ============================================================================
# RESULTS SUMMARY
# ============================================================================
cat("\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat("  ✓ ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

cat("📊 RESULTS SUMMARY:\n\n")
cat(sprintf("Cells analyzed:     %d\n", n_cells))
cat(sprintf("Genes detected:     %d\n", nrow(expr_matrix)))
cat(sprintf("Cell types:         %s\n", paste(unique(metadata$cell_type), collapse = ", ")))
cat(sprintf("Clusters found:     %d\n", n_distinct(metadata$cluster)))
cat(sprintf("Highly variable genes: %d\n", n_hvg))

cat("\n📁 FILES CREATED:\n\n")
cat("  ✓ data/processed/demo_analysis.rds - Complete analysis object\n")
cat("  ✓ data/processed/cell_metadata.csv - Cell annotations\n")
cat("  ✓ data/processed/marker_genes.csv - Top markers per cluster\n")
cat("  ✓ data/processed/pca_coordinates.csv - PCA coordinates\n")

cat("\n🚀 NEXT STEPS:\n\n")
cat("  1. Install Seurat for full scRNA-seq analysis:\n")
cat("     install.packages('Seurat')\n")
cat("     Then run: source('QUICKSTART.R')\n\n")
cat("  2. Launch Shiny dashboard:\n")
cat("     shiny::runApp('inst/app')\n\n")
cat("  3. Train ML models:\n")
cat("     source('ml_models.R')\n\n")
cat("  4. Integrate with LLM chatbot:\n")
cat("     source('chatbot_utils.R')\n\n")

cat("════════════════════════════════════════════════════════════════════════\n\n")
