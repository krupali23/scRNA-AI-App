# Data Processing Module for Single Cell RNA-Seq
# Handles loading, QC, and preprocessing of scRNA-seq data

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)

#' Load GSE Data Files
#' @param data_dir Directory containing raw scRNA-seq files
#' @param dataset_name Name of dataset (GSE123813 or GSE243013)
load_gse_data <- function(data_dir, dataset_name) {
  message(sprintf("Loading %s data...", dataset_name))
  
  if (dataset_name == "GSE123813") {
    # Load GSE123813 (BCC dataset)
    counts_file <- file.path(data_dir, "GSE123813_bcc_scRNA_counts.txt.gz")
    meta_file <- file.path(data_dir, "GSE123813_bcc_all_metadata.txt.gz")
    
    # Read counts matrix
    counts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
    
    # Read metadata
    metadata <- read.delim(meta_file, row.names = 1)
    
  } else if (dataset_name == "GSE243013") {
    # Load GSE243013 (NSCLC immune dataset)
    barcodes_file <- file.path(data_dir, "GSE243013_barcodes.csv.gz")
    genes_file <- file.path(data_dir, "GSE243013_genes.csv.gz")
    mtx_file <- file.path(data_dir, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
    meta_file <- file.path(data_dir, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
    
    # Read Matrix Market format
    library(Matrix)
    counts <- readMM(mtx_file) %>% as.matrix()
    
    genes <- read.csv(genes_file, header = FALSE, row.names = 1)
    barcodes <- read.csv(barcodes_file, header = FALSE)
    
    rownames(counts) <- rownames(genes)
    colnames(counts) <- barcodes$V1
    
    metadata <- read.csv(meta_file, row.names = 1)
  }
  
  message(sprintf("Loaded %d genes x %d cells", nrow(counts), ncol(counts)))
  return(list(counts = counts, metadata = metadata))
}

#' Create Seurat Object with QC
#' @param counts Gene expression matrix
#' @param metadata Cell metadata
create_seurat_object <- function(counts, metadata) {
  message("Creating Seurat object...")
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    min.cells = 3,
    min.features = 200
  )
  
  # QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )
  
  # Log normalize and scale
  seurat_obj <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData()
  
  # Dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
  
  # Clustering
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
  
  message(sprintf("Object created with %d cells, %d features",
                  ncol(seurat_obj), nrow(seurat_obj)))
  
  return(seurat_obj)
}

#' Perform Differential Expression Analysis
#' @param seurat_obj Processed Seurat object
#' @param ident1 First identity (cluster/celltype)
#' @param ident2 Second identity
perform_de <- function(seurat_obj, ident1, ident2) {
  message(sprintf("Running DE analysis: %s vs %s", ident1, ident2))
  
  de_results <- FindMarkers(
    seurat_obj,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  return(de_results %>% 
    rownames_to_column("gene") %>%
    arrange(p_val_adj))
}

#' Get cell type annotations with statistical summary
#' @param seurat_obj Processed Seurat object
get_celltype_summary <- function(seurat_obj) {
  df <- data.frame(
    cluster = Idents(seurat_obj),
    cell_id = colnames(seurat_obj),
    nUMI = colSums(seurat_obj@assays$RNA$counts),
    nGenes = colSums(seurat_obj@assays$RNA$counts > 0),
    percent_mt = seurat_obj$percent.mt
  )
  
  return(df)
}
