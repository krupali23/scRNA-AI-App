# Quick Start - Simplified (No optional dependencies)
# Works with just Seurat installed

cat("\n=== INSTALLING BIOCONDUCTOR PACKAGES ===\n")

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
bioc_packages <- c("SingleCellExperiment", "DESeq2", "scran", "scuttle")
for (pkg in bioc_packages) {
  if (!require(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("\n=== LOADING LIBRARIES ===\n")

library(Seurat)
library(tidyverse)
library(plotly)
library(DT)

cat("✓ All libraries loaded\n")

setwd("C:/Users/krupa/Desktop/Shiny_project")

# ============================================================================
# STEP 1: Load and Analyze Real Data
# ============================================================================

cat("\n=== STEP 1: LOADING REAL DATA ===\n")

data_dir <- "data/raw"

# Check what data files exist
data_files <- list.files(data_dir, pattern = "*.mtx|*.csv|*.tsv|*.txt")
cat(sprintf("Found %d data files\n", length(data_files)))
print(data_files)

# Load demo results if available (from previous run)
if (file.exists("data/processed/demo_analysis.rds")) {
  cat("\nLoading previous demo analysis results...\n")
  demo_data <- readRDS("data/processed/demo_analysis.rds")
  
  cat(sprintf("✓ Loaded demo analysis: %d cells, %d genes\n", 
              ncol(demo_data$matrix), nrow(demo_data$matrix)))
  
  # Create Seurat object from demo data
  if (is.data.frame(demo_data$matrix) || is.matrix(demo_data$matrix)) {
    seurat_obj <- CreateSeuratObject(counts = as.matrix(demo_data$matrix))
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 500)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 30)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
    
    cat("\n✓ Seurat object created and processed\n")
  }
}

# ============================================================================
# STEP 2: Clustering and Visualization
# ============================================================================

cat("\n=== STEP 2: CLUSTERING ===\n")

# Add clustering if object exists
if (exists("seurat_obj")) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  n_clusters <- length(unique(Idents(seurat_obj)))
  cat(sprintf("✓ Found %d clusters\n", n_clusters))
}

# ============================================================================
# STEP 3: ML Models
# ============================================================================

cat("\n=== STEP 3: ML MODELS ===\n")

if (exists("seurat_obj")) {
  # Simple classifier using top genes
  library(caret)
  library(randomForest)
  
  cat("Preparing ML training data...\n")
  
  # Get expression matrix
  expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
  
  # Get top variable genes
  top_genes <- VariableFeatures(seurat_obj)[1:100]
  expr_subset <- expr_matrix[top_genes, ]
  
  # Prepare data frame
  df_train <- data.frame(
    t(expr_subset),
    cluster = Idents(seurat_obj)
  )
  
  # Train RF model
  cat("Training Random Forest model (this may take a minute)...\n")
  set.seed(42)
  
  rf_model <- randomForest(
    cluster ~ .,
    data = df_train[sample(1:nrow(df_train), min(5000, nrow(df_train))), ],
    ntree = 100,
    importance = TRUE
  )
  
  accuracy <- mean(rf_model$predicted == df_train$cluster[
    which(rownames(df_train) %in% names(rf_model$predicted))
  ])
  
  cat(sprintf("✓ Random Forest model trained\n"))
  cat(sprintf("  OOB Error Rate: %.2f%%\n", rf_model$err.rate[nrow(rf_model$err.rate), 1] * 100))
  
  # Get feature importance
  importance_df <- data.frame(
    gene = rownames(rf_model$importance),
    importance = rf_model$importance[, 1]
  ) %>%
    arrange(desc(importance)) %>%
    head(20)
  
  cat("\nTop 20 Important Features:\n")
  print(importance_df)
  
  # Save model
  saveRDS(rf_model, "data/processed/ml_model_rf.rds")
  cat("✓ Saved model to data/processed/ml_model_rf.rds\n")
}

# ============================================================================
# STEP 4: GitHub Setup Check
# ============================================================================

cat("\n=== STEP 4: GITHUB SETUP ===\n")

if (file.exists("GITHUB_SETUP.md")) {
  cat("✓ GITHUB_SETUP.md exists\n")
  cat("  Follow steps in GITHUB_SETUP.md to publish to GitHub\n")
} else {
  cat("⚠️  GITHUB_SETUP.md not found\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n\n")

if (exists("seurat_obj")) {
  cat(sprintf("Seurat Object Summary:\n"))
  cat(sprintf("  Cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("  Genes: %d\n", nrow(seurat_obj)))
  cat(sprintf("  Clusters: %d\n", length(unique(Idents(seurat_obj)))))
  
  # Save for dashboard
  saveRDS(seurat_obj, "data/processed/seurat_full.rds")
  cat("\n✓ Saved full Seurat object to data/processed/seurat_full.rds\n")
}

cat("\nNext Steps:\n")
cat("  A) View results in dashboard: Open http://127.0.0.1:3838\n")
cat("  B) Publish to GitHub: Follow GITHUB_SETUP.md\n")
cat("  C) Explore full analysis in R: load('data/processed/seurat_full.rds')\n")

cat("\n")
