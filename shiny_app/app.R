library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(tidyverse)

# Load analysis results
analysis_data <- tryCatch({
  readRDS("../data/processed/demo_analysis.rds")
}, error = function(e) {
  NULL
})

marker_genes <- tryCatch({
  read.csv("../data/processed/marker_genes.csv")
}, error = function(e) {
  data.frame(cluster = character(), marker_gene = character())
})

pca_coords <- tryCatch({
  read.csv("../data/processed/pca_coordinates.csv")
}, error = function(e) {
  NULL
})

cell_metadata <- tryCatch({
  read.csv("../data/processed/cell_metadata.csv")
}, error = function(e) {
  NULL
})

# ============================================================================
# UI
# ============================================================================
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = "scRNA-seq Analysis Dashboard",
    tags$li(class = "dropdown", 
            HTML('<span style="padding:10px; color:white;">🧬 Bioinformatics</span>'))
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("home")),
      menuItem("Data Overview", tabName = "data", icon = icon("table")),
      menuItem("Quality Control", tabName = "qc", icon = icon("chart-bar")),
      menuItem("Dimensionality Reduction", tabName = "dimred", icon = icon("project-diagram")),
      menuItem("Clustering", tabName = "clustering", icon = icon("layer-group")),
      menuItem("Marker Genes", tabName = "markers", icon = icon("dna")),
      menuItem("ML Models", tabName = "ml", icon = icon("brain")),
      menuItem("AI Chatbot", tabName = "chat", icon = icon("comments")),
      menuItem("Documentation", tabName = "docs", icon = icon("book"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # ===== TAB 1: DASHBOARD =====
      tabItem(tabName = "dashboard",
        fluidRow(
          box(title = "Welcome", width = 12, status = "primary",
            h3("🧬 scRNA-seq Analysis Platform"),
            p("Interactive analysis dashboard with AI integration for single-cell RNA sequencing data."),
            br(),
            h4("Quick Stats:"),
            if (!is.null(cell_metadata)) {
              tags$ul(
                tags$li(sprintf("Total Cells: %d", nrow(cell_metadata))),
                tags$li(sprintf("Cell Types: %s", paste(unique(cell_metadata$cell_type), collapse = ", "))),
                tags$li(sprintf("Clusters: %d", length(unique(cell_metadata$cluster))))
              )
            } else {
              p("Loading data...")
            }
          )
        ),
        fluidRow(
          box(title = "Features", width = 6, status = "info",
            tags$ul(
              tags$li("Quality control metrics"),
              tags$li("Dimensionality reduction (PCA, UMAP, t-SNE)"),
              tags$li("Clustering analysis"),
              tags$li("Marker gene identification"),
              tags$li("ML model training"),
              tags$li("AI chatbot integration")
            )
          ),
          box(title = "Getting Started", width = 6, status = "success",
            tags$ol(
              tags$li("Explore raw data in 'Data Overview'"),
              tags$li("Check quality metrics in 'QC'"),
              tags$li("View clusters in 'Clustering'"),
              tags$li("Identify markers in 'Marker Genes'"),
              tags$li("Train ML models in 'ML Models'"),
              tags$li("Ask the AI chatbot questions")
            )
          )
        )
      ),
      
      # ===== TAB 2: DATA OVERVIEW =====
      tabItem(tabName = "data",
        fluidRow(
          box(title = "Cell Metadata", width = 12, status = "primary",
            DTOutput("metadata_table"),
            downloadButton("download_metadata", "Download CSV")
          )
        ),
        fluidRow(
          box(title = "Cell Type Distribution", width = 6,
            plotlyOutput("celltype_barplot")
          ),
          box(title = "Batch Distribution", width = 6,
            plotlyOutput("batch_barplot")
          )
        )
      ),
      
      # ===== TAB 3: QUALITY CONTROL =====
      tabItem(tabName = "qc",
        fluidRow(
          box(title = "Counts per Cell", width = 6,
            plotlyOutput("counts_histogram")
          ),
          box(title = "Genes per Cell", width = 6,
            plotlyOutput("genes_histogram")
          )
        ),
        fluidRow(
          box(title = "Counts vs Genes", width = 12,
            plotlyOutput("counts_vs_genes")
          )
        )
      ),
      
      # ===== TAB 4: DIMENSIONALITY REDUCTION =====
      tabItem(tabName = "dimred",
        fluidRow(
          box(title = "PCA Plot (Colored by Cell Type)", width = 12,
            plotlyOutput("pca_plot")
          )
        ),
        fluidRow(
          box(title = "PCA Plot (Colored by Cluster)", width = 12,
            plotlyOutput("pca_cluster_plot")
          )
        )
        ,
        fluidRow(
          box(title = "UMAP Previews (per-resolution)", width = 12, status = "primary",
              fluidRow(
                column(width = 3,
                    selectInput("umap_selector", "Select UMAP image:", choices = NULL),
                    downloadButton("download_umap", "Download PNG"),
                    br(),
                    tags$div(style='margin-top:8px;',
                        tags$strong('Gallery (click to select):'),
                        uiOutput('umap_gallery')
                    )
                ),
                column(width = 9,
                       uiOutput("umap_image_ui")
                )
              )
          )
        )
      ),
      
      # ===== TAB 5: CLUSTERING =====
      tabItem(tabName = "clustering",
        fluidRow(
          box(title = "Cluster Distribution", width = 6,
            plotlyOutput("cluster_barplot")
          ),
          box(title = "Cluster Statistics", width = 6,
            DTOutput("cluster_stats_table")
          )
        ),
        fluidRow(
          box(title = "Cell Type Composition by Cluster", width = 12,
            plotlyOutput("cluster_composition")
          )
        )
      ),
      
      # ===== TAB 6: MARKER GENES =====
      tabItem(tabName = "markers",
        fluidRow(
          box(title = "Marker Genes Table", width = 12, status = "primary",
            DTOutput("markers_table"),
            downloadButton("download_markers", "Download CSV")
          )
        ),
        fluidRow(
          box(title = "Top Markers per Cluster", width = 12,
            plotlyOutput("markers_barplot")
          )
        )
      ),
      
      # ===== TAB 7: ML MODELS =====
      tabItem(tabName = "ml",
        fluidRow(
          box(title = "Model Training", width = 12, status = "warning",
            p("Machine Learning Models (Coming Soon)"),
            tags$ul(
              tags$li("Random Forest Classifier"),
              tags$li("XGBoost Model"),
              tags$li("Feature Importance (SHAP)"),
              tags$li("Cross-validation Results")
            ),
            p("Once Seurat and caret packages are installed, you can train cell type prediction models."),
            actionButton("train_btn", "Train Models", class = "btn-primary", style = "margin-top: 10px;")
          )
        )
      ),
      
      # ===== TAB 8: AI CHATBOT =====
      tabItem(tabName = "chat",
        fluidRow(
          box(title = "AI Chatbot (LLM Integration)", width = 12, status = "info",
            p("Ask biological questions about your data:"),
            textInput("chat_input", "Your question:", placeholder = "e.g., What are the main cell types in this dataset?"),
            actionButton("chat_btn", "Send", class = "btn-primary"),
            hr(),
            verbatimTextOutput("chat_response")
          )
        )
      ),
      
      # ===== TAB 9: DOCUMENTATION =====
      tabItem(tabName = "docs",
        fluidRow(
          box(title = "Documentation", width = 12, status = "success",
            h4("Project Structure:"),
            tags$pre(
              "Shiny_project/
  ├── data/
  │   ├── raw/          # Raw sequencing files (7.2 GB)
  │   └── processed/    # Analysis results
  ├── ml_models/        # Trained models
  ├── shiny_app/        # This dashboard
  ├── data_processing.R # scRNA-seq pipeline
  ├── ml_models.R       # ML classification
  ├── chatbot_utils.R   # LLM integration
  ├── RUN_DEMO.R        # Demo analysis
  ├── README.md         # Full documentation
  └── GITHUB_SETUP.md   # Publishing guide"
            ),
            h4("Next Steps:"),
            tags$ol(
              tags$li("Install Seurat: install.packages('Seurat')"),
              tags$li("Run full analysis: source('QUICKSTART.R')"),
              tags$li("Configure LLM API keys for chatbot"),
              tags$li("Publish to GitHub for portfolio")
            ),
            h4("Files Created:"),
            tags$ul(
              tags$li("data/processed/demo_analysis.rds - Analysis object"),
              tags$li("data/processed/cell_metadata.csv - Cell annotations"),
              tags$li("data/processed/marker_genes.csv - Top markers"),
              tags$li("data/processed/pca_coordinates.csv - PCA data")
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {
  # Expose per-resolution UMAP images directory to the Shiny static resources
  umap_dir <- normalizePath(file.path("..","data","processed","umap_per_resolution"), mustWork = FALSE)
  if (dir.exists(umap_dir)) {
    addResourcePath('umap_images', umap_dir)
  }

  # Helper to list available UMAP images (basename only)
  list_umap_images <- function(){
    if(dir.exists(umap_dir)){
      tools::file_path_sans_ext(basename(list.files(umap_dir, pattern='\\.png$', full.names=TRUE)))
    } else character(0)
  }

  # Populate selector once on app start
  observe({
    imgs <- list.files(umap_dir, pattern='\\.png$', full.names=FALSE)
    if(length(imgs)>0){
      updateSelectInput(session, 'umap_selector', choices = imgs, selected = imgs[1])
    }
  })

  # Render UMAP image UI (uses the resource path 'umap_images')
  output$umap_image_ui <- renderUI({
    imgs <- list.files(umap_dir, pattern='\\.png$', full.names=FALSE)
    if(length(imgs)==0) return(tags$p('No UMAP images found in data/processed/umap_per_resolution'))
    sel <- input$umap_selector
    if(is.null(sel) || !(sel %in% imgs)) sel <- imgs[1]
    tags$div(style = 'text-align:center;',
             tags$img(src = file.path('umap_images', sel), style = 'max-width:100%; height:auto; border:1px solid #ddd;'))
  })

  # Render clickable gallery of thumbnails
  output$umap_gallery <- renderUI({
    imgs <- list.files(umap_dir, pattern='\\.png$', full.names=FALSE)
    if(length(imgs)==0) return(NULL)
    # create clickable thumbnails that set the select input when clicked
    thumbs <- lapply(imgs, function(fn){
      # safe JS to set input value
      js <- sprintf("Shiny.setInputValue('umap_selector', '%s', {priority: 'event'})", fn)
      tags$div(style='display:inline-block; margin:4px; text-align:center;',
               tags$img(src = file.path('umap_images', fn), onclick = HTML(js), style='width:120px; height:auto; cursor:pointer; border:1px solid #ccc;'),
               tags$div(style='font-size:11px; margin-top:4px; max-width:120px; word-wrap:break-word;', fn)
      )
    })
    tags$div(style='display:flex; flex-wrap:wrap; align-items:flex-start;', thumbs)
  })

  # Download selected UMAP image
  output$download_umap <- downloadHandler(
    filename = function(){
      if(is.null(input$umap_selector) || input$umap_selector=='') return('umap.png')
      input$umap_selector
    },
    content = function(file){
      src <- file.path(umap_dir, input$umap_selector)
      file.copy(src, file, overwrite=TRUE)
    }
  )
  
  # Reactive data
  metadata_reactive <- reactive({
    if (!is.null(cell_metadata)) {
      head(cell_metadata, 100)
    } else {
      data.frame()
    }
  })
  
  # TABLE: Metadata
  output$metadata_table <- renderDT({
    datatable(metadata_reactive(), options = list(pageLength = 10))
  })
  
  # PLOT: Cell type barplot
  output$celltype_barplot <- renderPlotly({
    if (!is.null(cell_metadata)) {
      df <- cell_metadata %>%
        group_by(cell_type) %>%
        summarise(count = n(), .groups = 'drop') %>%
        arrange(desc(count))
      
      plot_ly(df, x = ~cell_type, y = ~count, type = 'bar', marker = list(color = '#1f77b4')) %>%
        layout(title = "Cell Type Distribution", xaxis = list(title = "Cell Type"), yaxis = list(title = "Count"))
    }
  })
  
  # PLOT: Batch barplot
  output$batch_barplot <- renderPlotly({
    if (!is.null(cell_metadata)) {
      df <- cell_metadata %>%
        group_by(batch) %>%
        summarise(count = n(), .groups = 'drop')
      
      plot_ly(df, x = ~batch, y = ~count, type = 'bar', marker = list(color = '#ff7f0e')) %>%
        layout(title = "Batch Distribution", xaxis = list(title = "Batch"), yaxis = list(title = "Count"))
    }
  })
  
  # PLOT: Counts histogram
  output$counts_histogram <- renderPlotly({
    if (!is.null(cell_metadata)) {
      plot_ly(x = cell_metadata$n_counts, type = 'histogram', nbinsx = 30, marker = list(color = '#2ca02c')) %>%
        layout(title = "Distribution of Counts per Cell", xaxis = list(title = "n_counts"), yaxis = list(title = "Frequency"))
    }
  })
  
  # PLOT: Genes histogram
  output$genes_histogram <- renderPlotly({
    if (!is.null(cell_metadata)) {
      plot_ly(x = cell_metadata$n_genes, type = 'histogram', nbinsx = 30, marker = list(color = '#d62728')) %>%
        layout(title = "Distribution of Genes per Cell", xaxis = list(title = "n_genes"), yaxis = list(title = "Frequency"))
    }
  })
  
  # PLOT: Counts vs Genes
  output$counts_vs_genes <- renderPlotly({
    if (!is.null(cell_metadata)) {
      plot_ly(cell_metadata, x = ~n_genes, y = ~n_counts, 
              mode = 'markers', type = 'scatter', marker = list(size = 5, color = ~n_counts, colorscale = 'Viridis'),
              text = ~cell_type, hoverinfo = 'text') %>%
        layout(title = "Genes vs Counts per Cell", xaxis = list(title = "n_genes"), yaxis = list(title = "n_counts"))
    }
  })
  
  # PLOT: PCA (Cell Type)
  output$pca_plot <- renderPlotly({
    if (!is.null(pca_coords)) {
      plot_ly(pca_coords, x = ~PC1, y = ~PC2, color = ~cell_type, type = 'scatter', mode = 'markers',
              marker = list(size = 8), colors = "Set1") %>%
        layout(title = "PCA Plot - Colored by Cell Type", xaxis = list(title = "PC1"), yaxis = list(title = "PC2"))
    }
  })
  
  # PLOT: PCA (Cluster)
  output$pca_cluster_plot <- renderPlotly({
    if (!is.null(pca_coords)) {
      plot_ly(pca_coords, x = ~PC1, y = ~PC2, color = ~cluster, type = 'scatter', mode = 'markers',
              marker = list(size = 8), colors = "Set2") %>%
        layout(title = "PCA Plot - Colored by Cluster", xaxis = list(title = "PC1"), yaxis = list(title = "PC2"))
    }
  })
  
  # TABLE: Cluster stats
  output$cluster_stats_table <- renderDT({
    if (!is.null(cell_metadata)) {
      stats <- cell_metadata %>%
        group_by(cluster) %>%
        summarise(
          n_cells = n(),
          dominant_type = names(which.max(table(cell_type))),
          n_types = n_distinct(cell_type),
          mean_counts = round(mean(n_counts), 1),
          .groups = 'drop'
        )
      datatable(stats, options = list(pageLength = 10))
    }
  })
  
  # PLOT: Cluster barplot
  output$cluster_barplot <- renderPlotly({
    if (!is.null(cell_metadata)) {
      df <- cell_metadata %>%
        group_by(cluster) %>%
        summarise(count = n(), .groups = 'drop') %>%
        arrange(as.numeric(cluster))
      
      plot_ly(df, x = ~cluster, y = ~count, type = 'bar', marker = list(color = '#9467bd')) %>%
        layout(title = "Cell Distribution by Cluster", xaxis = list(title = "Cluster"), yaxis = list(title = "Count"))
    }
  })
  
  # PLOT: Cluster composition
  output$cluster_composition <- renderPlotly({
    if (!is.null(cell_metadata)) {
      df <- cell_metadata %>%
        group_by(cluster, cell_type) %>%
        summarise(count = n(), .groups = 'drop')
      
      plot_ly(df, x = ~cluster, y = ~count, fill = ~cell_type, type = 'bar') %>%
        layout(barmode = 'stack', title = "Cell Type Composition by Cluster", 
               xaxis = list(title = "Cluster"), yaxis = list(title = "Count"))
    }
  })
  
  # TABLE: Marker genes
  output$markers_table <- renderDT({
    datatable(marker_genes, options = list(pageLength = 15))
  })
  
  # PLOT: Top markers
  output$markers_barplot <- renderPlotly({
    if (nrow(marker_genes) > 0) {
      top_markers <- marker_genes %>%
        group_by(cluster) %>%
        slice_head(n = 5) %>%
        ungroup()
      
      plot_ly(top_markers, x = ~marker_gene, y = ~cluster, type = 'bar', orientation = 'h', marker = list(color = '#17becf')) %>%
        layout(title = "Top 5 Markers per Cluster", xaxis = list(title = "Marker Gene"), yaxis = list(title = "Cluster"))
    }
  })
  
  # DOWNLOAD: Metadata
  output$download_metadata <- downloadHandler(
    filename = "cell_metadata.csv",
    content = function(file) {
      write.csv(cell_metadata, file, row.names = FALSE)
    }
  )
  
  # DOWNLOAD: Markers
  output$download_markers <- downloadHandler(
    filename = "marker_genes.csv",
    content = function(file) {
      write.csv(marker_genes, file, row.names = FALSE)
    }
  )
  
  # CHATBOT
  observeEvent(input$chat_btn, {
    question <- input$chat_input
    
    # Simple response logic
    response <- if (grepl("cell type", question, ignore.case = TRUE)) {
      sprintf("This dataset contains %d cells from %d cell types: %s. The dominant cell type is %s.",
              nrow(cell_metadata),
              n_distinct(cell_metadata$cell_type),
              paste(unique(cell_metadata$cell_type), collapse = ", "),
              names(which.max(table(cell_metadata$cell_type))))
    } else if (grepl("cluster", question, ignore.case = TRUE)) {
      sprintf("The analysis identified %d distinct clusters using k-means clustering on the top 10 PCs.",
              n_distinct(cell_metadata$cluster))
    } else if (grepl("marker", question, ignore.case = TRUE)) {
      sprintf("Top marker genes were identified using log fold change comparison. See the Marker Genes tab for details.")
    } else {
      "I'm a simple chatbot for this dataset. Ask about cell types, clusters, or markers to get started!"
    }
    
    output$chat_response <- renderText({
      paste("🤖 Assistant:\n", response)
    })
  })
  
  # ML Models button
  observeEvent(input$train_btn, {
    showNotification("Model training requires Seurat and caret packages. Install with: install.packages(c('Seurat', 'caret', 'randomForest', 'xgboost'))",
                     type = "info", duration = 5)
  })
}

# ============================================================================
# RUN APP
# ============================================================================
shinyApp(ui, server)
