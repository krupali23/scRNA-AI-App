#!/usr/bin/env Rscript
options(repos='https://cloud.r-project.org')
.libPaths(c('C:/Users/krupa/AppData/Local/R/library', .libPaths()))

safe_load <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      install.packages(p, dependencies=TRUE)
    }
    library(p, character.only=TRUE)
  }
}

safe_load(c('Seurat','ggplot2','readr','dplyr','patchwork'))

out_dir <- file.path('data','processed')
dir.create(file.path(out_dir,'umap_per_resolution'), recursive=TRUE, showWarnings=FALSE)

message('Loading Seurat object (mapped preferred)...')
seurat_path <- file.path(out_dir,'seurat_full_mapped.rds')
if(!file.exists(seurat_path)) seurat_path <- file.path(out_dir,'seurat_full.rds')
if(!file.exists(seurat_path)) stop('No Seurat object found in data/processed/')
seu <- readRDS(seurat_path)
message('Loaded: ', seurat_path)

if(!'umap' %in% names(seu@reductions)){
  if('pca' %in% names(seu@reductions)){
    message('UMAP not found; attempting RunUMAP using first 20 PCs...')
    seu <- tryCatch({
      seu <- RunUMAP(seu, reduction='pca', dims=1:20)
    }, error=function(e){
      message('RunUMAP failed: ', e$message)
      return(seu)
    })
  }
}

# Find cluster files
cl_files <- list.files(out_dir, pattern='clusters_resolution_.*\\.csv$', full.names=TRUE)
if(length(cl_files)==0){
  message('No per-resolution cluster CSVs found. Looking for cluster summary...')
  sumf <- file.path(out_dir,'cluster_summary_resolutions.csv')
  if(file.exists(sumf)){
    df <- readr::read_csv(sumf, show_col_types=FALSE)
    if('resolution' %in% names(df)){
      res <- df$resolution
      # try to find matching csv files by resolution
      for(r in res){
        pattern <- paste0('clusters_resolution_', r)
        matches <- list.files(out_dir, pattern=pattern, full.names=TRUE)
        cl_files <- c(cl_files, matches)
      }
    }
  }
}

plots <- list()
if(length(cl_files)==0){
  message('No cluster files found. Using existing identity classes if present.')
  if(!is.null(Idents(seu))){
    p <- DimPlot(seu, reduction='umap', group.by=NULL) + ggtitle('UMAP (Idents)')
    ggsave(file.path(out_dir,'umap_per_resolution','umap_Idents.png'), p, width=6, height=5, dpi=150)
    message('Saved: ', file.path(out_dir,'umap_per_resolution','umap_Idents.png'))
  }
} else {
  for(f in cl_files){
    message('Processing clusters file: ', f)
    cf <- readr::read_csv(f, show_col_types=FALSE)
    # heuristics: find column with cells and cluster
    cn <- names(cf)
    # common formats: cell, cluster OR barcode, cluster
    cell_col <- cn[grepl('cell|barcode|cell_barcode', cn, ignore.case=TRUE)][1]
    cluster_col <- cn[grepl('cluster|assignment|cluster_id', cn, ignore.case=TRUE)][1]
    if(is.na(cell_col) || is.na(cluster_col)){
      # fallback: first two columns
      cell_col <- cn[1]; cluster_col <- cn[2]
    }
    # ensure matching names as characters
    cf[[cell_col]] <- as.character(cf[[cell_col]])
    cf[[cluster_col]] <- as.character(cf[[cluster_col]])
    # match to Seurat columns
    # seurat cellnames are colnames(seu)
    meta <- seu@meta.data
    meta$.__cellname__ <- rownames(meta)
    # try join by cell barcode
    joined <- tryCatch({
      left_join(data.frame(.__cellname__=colnames(seu)), cf, by=setNames(cell_col, '.__cellname__'))
    }, error=function(e){
      # try different approach: if cf has rownames
      if(nrow(cf)==ncol(seu)){
        df2 <- cf
        df2$.__cellname__ <- colnames(seu)
        df2
      } else stop('Could not align cluster file to Seurat object')
    })
    if(!('.__cellname__' %in% names(joined))) stop('Alignment failed')
    colname_meta <- paste0('cluster_res_', tools::file_path_sans_ext(basename(f)))
    # use cluster_col values
    # if joined has the cluster column name, extract it
    if(cluster_col %in% names(joined)){
      newcol <- joined[[cluster_col]]
    } else {
      newcol <- joined[[2]]
    }
    # assign to Seurat meta
    seu@meta.data[[colname_meta]] <- as.factor(newcol)
    p <- DimPlot(seu, group.by=colname_meta, reduction='umap') + ggtitle(colname_meta)
    outfn <- file.path(out_dir,'umap_per_resolution', paste0(colname_meta, '.png'))
    ggsave(outfn, p, width=6, height=5, dpi=150)
    message('Saved: ', outfn)
    plots[[colname_meta]] <- p
  }
  # save combined
  if(length(plots)>0){
    comb <- wrap_plots(plots, ncol=2)
    ggsave(file.path(out_dir,'umap_per_resolution','umap_per_resolutions_combined.png'), comb, width=10, height=5*ceiling(length(plots)/2), dpi=150)
    message('Saved combined UMAP: ', file.path(out_dir,'umap_per_resolution','umap_per_resolutions_combined.png'))
  }
}

# ---- Marker overlap analysis ----
message('\nRunning marker overlap analysis...')
top20f <- file.path(out_dir,'top20_markers_per_cluster_mapped.csv')
if(!file.exists(top20f)) top20f <- file.path(out_dir,'top20_markers_per_cluster.csv')
if(!file.exists(top20f)) stop('No top20 markers CSV found in data/processed/')
top20 <- readr::read_csv(top20f, show_col_types=FALSE)

# collect published gene lists from data/raw
pub_files <- list.files('data/raw', pattern='(marker|gene[_-]?list|published|signature).*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if(length(pub_files)==0){
  message('No published gene list CSVs found under data/raw. Will try to use demo markers if available.')
}

overlaps <- list()
if(length(pub_files)>0){
  for(pf in pub_files){
    message('Comparing to published list: ', pf)
    pub <- tryCatch({
      readr::read_csv(pf, show_col_types=FALSE)
    }, error=function(e){
      readr::read_delim(pf, delim='\\t', show_col_types=FALSE)
    })
    # find first column with gene symbols
    gcol <- names(pub)[which.max(vapply(pub, function(x) sum(!is.na(x)), integer(1)))]
    pub_genes <- unique(as.character(pub[[gcol]]))
    pub_genes <- pub_genes[!is.na(pub_genes) & pub_genes!='']
    # for each cluster
    clusters <- unique(top20$cluster)
    for(cl in clusters){
      tx <- top20 %>% filter(cluster==cl) %>% pull(gene)
      ov <- intersect(tx, pub_genes)
      overlaps[[length(overlaps)+1]] <- data.frame(published_file=basename(pf), cluster=cl, n_top20=length(tx), n_overlap=length(ov), pct_overlap=round(100*length(ov)/max(1,length(tx)),2), overlapping_genes=paste(ov, collapse=';'))
    }
  }
} else {
  # fallback: compare to demo_data markers if available inside data/processed/demo_analysis.rds
  demo_rds <- file.path(out_dir,'demo_analysis.rds')
  if(file.exists(demo_rds)){
    demo <- readRDS(demo_rds)
    if('marker_genes' %in% names(demo)){
      pub_genes <- unlist(demo$marker_genes)
      clusters <- unique(top20$cluster)
      for(cl in clusters){
        tx <- top20 %>% filter(cluster==cl) %>% pull(gene)
        ov <- intersect(tx, pub_genes)
        overlaps[[length(overlaps)+1]] <- data.frame(published_file='demo_marker_list', cluster=cl, n_top20=length(tx), n_overlap=length(ov), pct_overlap=round(100*length(ov)/max(1,length(tx)),2), overlapping_genes=paste(ov, collapse=';'))
      }
    } else message('demo_analysis.rds has no marker_genes element; skipping fallback overlap')
  } else message('No demo_analysis.rds for fallback; skipping overlap analysis')
}

if(length(overlaps)>0){
  ovdf <- do.call(rbind, overlaps)
  out_csv <- file.path(out_dir,'marker_overlap_summary.csv')
  readr::write_csv(ovdf, out_csv)
  message('Wrote overlap summary to ', out_csv)
  # short text report
  rpt <- file.path(out_dir,'marker_overlap_report.txt')
  sink(rpt)
  cat('Marker overlap report\n')
  cat('====================\n\n')
  cat('Seurat object: ', seurat_path, '\n')
  cat('Top20 markers source: ', top20f, '\n\n')
  for(pf in unique(ovdf$published_file)){
    sub <- ovdf %>% filter(published_file==pf)
    cat('Published set: ', pf, '\n')
    cat('Clusters with highest overlap:\n')
    disp <- sub %>% arrange(desc(n_overlap)) %>% head(5)
    print(disp)
    cat('\n')
  }
  sink()
  message('Wrote textual report to ', rpt)
} else {
  message('No overlaps computed.')
}

message('\nAll done. UMAPs and overlap files saved under data/processed/.')
