# Quick Start - Simplified (No optional dependencies)
# Works with just Seurat installed

cat("\n=== INSTALLING BIOCONDUCTOR PACKAGES ===\n")

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
## On Windows we skip source-only packages (e.g. 'bluster') to avoid
## compilation failures. Only install essential Bioconductor packages here.
bioc_packages <- c("SingleCellExperiment", "DESeq2")
for (pkg in bioc_packages) {
  if (!require(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE),
      error = function(e) {
        cat(sprintf("Warning: installation of %s failed: %s\n", pkg, e$message))
      }
    )
  }
}

cat("Note: 'scran' and 'scuttle' installation skipped in simplified run to avoid Windows-only compilation issues.\n")

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

  # If the saved object is already a Seurat object, reuse it
  if (inherits(demo_data, "Seurat")) {
    seurat_obj <- demo_data
    cat("✓ demo_analysis.rds contains a Seurat object; reusing it.\n")
  } else if (inherits(demo_data, "SingleCellExperiment") || inherits(demo_data, "SummarizedExperiment")) {
    # Try to extract counts/assay
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      counts_mat <- tryCatch(SummarizedExperiment::assay(demo_data), error = function(e) NULL)
    } else {
      counts_mat <- NULL
    }
    if (!is.null(counts_mat)) {
      cat(sprintf("✓ Extracted assay matrix from SingleCellExperiment/SummarizedExperiment: %d cells, %d genes\n", ncol(counts_mat), nrow(counts_mat)))
      seurat_obj <- CreateSeuratObject(counts = as.matrix(counts_mat))
    }
  } else {
    # Look for common fields that may hold the expression/counts matrix
    candidate_names <- c("expression_matrix", "normalized_matrix", "matrix", "counts", "expr", "expression", "data")
    found <- FALSE
    for (nm in candidate_names) {
      if (!is.null(demo_data[[nm]])) {
        mat <- demo_data[[nm]]
        if (is.data.frame(mat) || is.matrix(mat)) {
          cat(sprintf("✓ Found demo_data[['%s']] with dimensions %dx%d\n", nm, ncol(mat), nrow(mat)))
          seurat_obj <- CreateSeuratObject(counts = as.matrix(mat))
          found <- TRUE
          break
        }
      }
    }
    if (!found) {
      cat("⚠️  demo_analysis.rds loaded but no suitable count matrix found (matrix/counts/expr).\n")
    }
  }

  if (exists("seurat_obj")) {
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

# Clustering parameters (tweak these for more/fewer clusters)
cluster_resolution <- 1.0
neighbor_dims <- 1:30

# Add clustering if object exists
if (exists("seurat_obj")) {
  cat(sprintf("Running FindNeighbors(dims = %s) and FindClusters(res = %.2f)\n",
              paste(range(neighbor_dims), collapse = ":"), cluster_resolution))
  seurat_obj <- FindNeighbors(seurat_obj, dims = neighbor_dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution)

  # Report cluster sizes
  cl_tab <- table(Idents(seurat_obj))
  cat("Cluster sizes:\n")
  print(cl_tab)

  n_clusters <- length(unique(Idents(seurat_obj)))
  cat(sprintf("✓ Found %d clusters\n", n_clusters))
}

# ============================================================================
# Multi-resolution clustering: run several resolutions and save assignments
# ============================================================================
if (exists("seurat_obj")) {
  cat("\n=== MULTI-RESOLUTION CLUSTERING ===\n")
  resolutions <- c(0.2, 0.5, 0.8, 1.0, 1.5)
  cluster_summary <- data.frame(resolution = numeric(), n_clusters = integer(), stringsAsFactors = FALSE)

  for (res in resolutions) {
    cat(sprintf("Running clustering at resolution = %.2f\n", res))
    tryCatch({
      seurat_obj <- FindClusters(seurat_obj, resolution = res)
      col_name <- paste0('cluster_res_', gsub("\\.", "-", as.character(res)))
      # store cluster assignments in meta.data under a stable column name
      seurat_obj@meta.data[[col_name]] <- as.character(Idents(seurat_obj))

      # write per-resolution cluster assignments
      out_df <- data.frame(cell = rownames(seurat_obj@meta.data), cluster = seurat_obj@meta.data[[col_name]], stringsAsFactors = FALSE)
      out_path <- file.path('data', 'processed', sprintf('clusters_resolution_%s.csv', gsub('\\.', '-', as.character(res))))
      write.csv(out_df, out_path, row.names = FALSE)
      cat(sprintf("Saved cluster assignments to %s\n", out_path))

      # record summary
      ncl <- length(unique(out_df$cluster))
      cluster_summary <- rbind(cluster_summary, data.frame(resolution = res, n_clusters = ncl))
    }, error = function(e) {
      cat(sprintf("Warning: clustering failed at resolution %.2f: %s\n", res, e$message))
    })
  }

  # Save summary table
  summary_path <- file.path('data', 'processed', 'cluster_summary_resolutions.csv')
  write.csv(cluster_summary, summary_path, row.names = FALSE)
  cat(sprintf("Saved cluster summary to %s\n", summary_path))
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

  # Determine labels: prefer Idents(seurat_obj), fallback to demo_data clusters if available
  labels <- NULL
  if (length(unique(Idents(seurat_obj))) > 1) {
    labels <- Idents(seurat_obj)
  } else if (exists("demo_data") && !is.null(demo_data$clusters)) {
    labels <- demo_data$clusters
    cat("Using labels from demo_data$clusters for ML training.\n")
  } else if (exists("demo_data") && !is.null(demo_data$cluster)) {
    labels <- demo_data$cluster
    cat("Using labels from demo_data$cluster for ML training.\n")
  }

  if (is.null(labels) || length(unique(labels)) < 2) {
    cat("⚠️  Skipping ML: need at least two classes for classification.\n")
  } else {
    # Get top variable genes
    top_genes <- VariableFeatures(seurat_obj)[1:100]
    expr_subset <- expr_matrix[top_genes, ]

    # Prepare data frame
    df_train <- data.frame(
      t(expr_subset),
      cluster = as.factor(labels)
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
