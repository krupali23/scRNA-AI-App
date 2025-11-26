# Quick Start Example
# Load and analyze scRNA-seq data with AI chatbot

# Set working directory
setwd("C:/Users/krupa/Desktop/Shiny_project")

# Load required libraries
library(Seurat)
library(tidyverse)
library(plotly)
library(DT)

# Source utility functions (non-blocking: optional modules)
tryCatch(
  source("data_processing.R"),
  error = function(e) message("Note: data_processing.R failed to load (optional): ", e$message)
)
tryCatch(
  source("chatbot_utils.R"),
  error = function(e) message("Note: chatbot_utils.R failed to load (optional): ", e$message)
)
tryCatch(
  source("ml_models.R"),
  error = function(e) message("Note: ml_models.R failed to load (optional): ", e$message)
)

# ============================================================================
# STEP 1: Load Data
# ============================================================================

message("\n=== STEP 1: LOADING DATA ===\n")

data_dir <- "data/raw"

# Load GSE243013 (NSCLC Immune)
gse_data <- load_gse_data(data_dir, "GSE243013")

message("✓ Data loaded successfully!")
message(sprintf("  Genes: %d", nrow(gse_data$counts)))
message(sprintf("  Cells: %d", ncol(gse_data$counts)))

# ============================================================================
# STEP 2: Create Seurat Object with QC
# ============================================================================

message("\n=== STEP 2: CREATING SEURAT OBJECT ===\n")

seurat_obj <- create_seurat_object(gse_data$counts, gse_data$metadata)

message("✓ Seurat object created!")
message(sprintf("  Dimensions: %d x %d", nrow(seurat_obj), ncol(seurat_obj)))
message(sprintf("  Clusters: %d", length(unique(Idents(seurat_obj)))))

# ============================================================================
# STEP 3: Visualize Results
# ============================================================================

message("\n=== STEP 3: VISUALIZATION ===\n")

# UMAP plot
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  theme(legend.position = "right")

# Feature plot (CD4, CD8, CD14 - immune markers)
feature_plot <- FeaturePlot(seurat_obj, features = c("CD4", "CD8", "CD14"), 
                            reduction = "umap")

# QC metrics plot
qc_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3)

message("✓ Plots generated!")

# ============================================================================
# STEP 4: Find Markers & Analyze
# ============================================================================

message("\n=== STEP 4: FINDING MARKERS ===\n")

all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  arrange(cluster, desc(avg_log2FC))

message("✓ Top 10 markers per cluster:")
print(top_markers %>% select(cluster, gene, avg_log2FC, p_val_adj))

# ============================================================================
# STEP 5: Initialize AI Chatbot (Optional)
# ============================================================================

message("\n=== STEP 5: AI CHATBOT SETUP ===\n")

# Configure LLM (if API keys available)
if (Sys.getenv("AZURE_API_KEY") != "") {
  llm_config <- initialize_llm(
    provider = "azure",
    endpoint = Sys.getenv("AZURE_ENDPOINT"),
    api_key = Sys.getenv("AZURE_API_KEY")
  )
  
  message("✓ LLM configured (Azure OpenAI)")
  
  # Create RAG context
  rag_context <- create_rag_context(seurat_obj, 
                                     query_genes = c("CD4", "CD8", "CD14"))
  
  message("✓ RAG context created")
  
  # Example query
  query <- "What are the main immune cell types in this dataset?"
  # response <- query_llm(query, rag_context, llm_config)
  # message(sprintf("\nLLM Response: %s", response$response))
  
} else {
  message("⚠️  AZURE_API_KEY not set - skipping LLM initialization")
  message("   Set environment variables to enable AI chatbot")
}

# ============================================================================
# STEP 6: Train ML Model (Optional)
# ============================================================================

message("\n=== STEP 6: ML MODEL TRAINING ===\n")

if (FALSE) {  # Set to TRUE to run (takes time)
  trained_model <- train_celltype_classifier(seurat_obj, test_split = 0.2)
  
  message(sprintf("✓ Model trained with accuracy: %.3f", 
                  trained_model$performance$rf_accuracy))
  
  # Get feature importance
  importance <- get_feature_importance(trained_model, method = "rf", top_n = 20)
  message("Top 20 important genes:")
  print(importance)
}

# ============================================================================
# STEP 7: Summary Statistics
# ============================================================================

message("\n=== ANALYSIS SUMMARY ===\n")

celltype_summary <- get_celltype_summary(seurat_obj)

summary_table <- celltype_summary %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_genes = mean(nGenes),
    mean_umi = mean(nUMI),
    mean_mt_pct = mean(percent_mt),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(cluster))

message("Cluster Summary:")
print(summary_table)

# ============================================================================
# SAVE RESULTS
# ============================================================================

message("\n=== SAVING RESULTS ===\n")

# Save Seurat object
saveRDS(seurat_obj, "data/processed/seurat_object.rds")
message("✓ Saved: data/processed/seurat_object.rds")

# Save markers
write.csv(top_markers, "data/processed/marker_genes.csv", row.names = FALSE)
message("✓ Saved: data/processed/marker_genes.csv")

# Save summary
write.csv(summary_table, "data/processed/cluster_summary.csv", row.names = FALSE)
message("✓ Saved: data/processed/cluster_summary.csv")

message("\n=== COMPLETE ===\n")
message("Next: Run the Shiny app with shiny::runApp('inst/app')")

# ============================================================================
# EXAMPLE: INTERACTIVE EXPLORATION (Uncomment to run in R console)
# ============================================================================

if (FALSE) {
  # View gene expression
  genes_of_interest <- c("CD4", "CD8A", "CD14", "FCER1G", "CD19", "FCGR3A")
  FeaturePlot(seurat_obj, features = genes_of_interest, ncol = 3)
  
  # Compare two clusters
  de_markers <- perform_de(seurat_obj, ident1 = "0", ident2 = "1")
  head(de_markers, 20)
  
  # Volcano plot
  plot(de_markers$avg_log2FC, -log10(de_markers$p_val_adj))
}
