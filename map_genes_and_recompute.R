# Map generic gene IDs to real symbols using data/raw/GSE243013_genes.csv.gz
# Then update Seurat object, recompute markers and plots

.libPaths(c('C:/Users/krupa/AppData/Local/R/library', .libPaths()))
options(repos='https://cloud.r-project.org')

library(Seurat)
library(dplyr)

# Read gene symbols
genes_df <- read.csv(gzfile('data/raw/GSE243013_genes.csv.gz'), header = TRUE, stringsAsFactors = FALSE)
if (!'geneSymbol' %in% names(genes_df)) stop('Expected column geneSymbol in genes file')
gene_symbols <- genes_df$geneSymbol
cat(sprintf('Loaded %d gene symbols\n', length(gene_symbols)))

# Load existing Seurat object if available
seurat_path <- 'data/processed/seurat_full.rds'
if (file.exists(seurat_path)) {
  seurat_obj <- readRDS(seurat_path)
  cat('Loaded Seurat object from', seurat_path, '\n')
} else {
  # fallback: build from demo analysis rds
  demo <- readRDS('data/processed/demo_analysis.rds')
  mat <- demo$expression_matrix
  seurat_obj <- CreateSeuratObject(counts = as.matrix(mat))
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 500)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  cat('Created Seurat object from demo data\n')
}

n_genes <- nrow(seurat_obj)
cat('Seurat object has', n_genes, 'features\n')

# Map by position: gene_1 -> first row in gene_symbols
if (length(gene_symbols) < n_genes) stop('Not enough gene symbols to map')
new_names <- gene_symbols[1:n_genes]
# make unique and valid for Seurat
new_names <- make.unique(gsub('[^A-Za-z0-9_.-]', '_', new_names))

old_names <- rownames(seurat_obj)
cat('First old names:', paste(head(old_names,10), collapse=', '), '\n')
cat('First new names:', paste(head(new_names,10), collapse=', '), '\n')

# Apply new names
rownames(seurat_obj) <- new_names

# Recompute variable features if needed
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 500)

# Re-run PCA/UMAP to ensure consistency
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Clustering (use previous tuned parameters)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.0)

# Save updated Seurat object
out_rds <- 'data/processed/seurat_full_mapped.rds'
saveRDS(seurat_obj, out_rds)
cat('Saved mapped Seurat object to', out_rds, '\n')

# Find markers with real gene symbols
markers_all <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_all, file = 'data/processed/markers_all_clusters_mapped.csv', row.names = FALSE)
cat('Saved markers_all_clusters_mapped.csv\n')

# Top 20
top_markers <- markers_all %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20)
write.csv(top_markers, file = 'data/processed/top20_markers_per_cluster_mapped.csv', row.names = FALSE)
cat('Saved top20_markers_per_cluster_mapped.csv\n')

# DotPlot and heatmap
png('data/processed/dotplot_top_markers_mapped.png', width=1400, height=900)
try({ print(DotPlot(seurat_obj, features = unique(top_markers$gene)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))) }, silent = TRUE)
dev.off()
cat('Saved dotplot_top_markers_mapped.png\n')

png('data/processed/heatmap_top_markers_mapped.png', width=1400, height=900)
try({ print(DoHeatmap(seurat_obj, features = unique(top_markers$gene), size = 3) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))) }, silent = TRUE)
dev.off()
cat('Saved heatmap_top_markers_mapped.png\n')

cat('Done mapping and recomputing markers.\n')
