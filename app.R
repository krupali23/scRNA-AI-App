# Main Shiny/Golem App - Single Cell RNA-Seq Analysis with AI
# Author: Your Name
# Description: Interactive app for scRNA-seq analysis with integrated AI chatbot

library(shiny)
library(golem)
library(DT)
library(ggplot2)
library(plotly)
library(shinydashboard)
library(shinyalert)
library(tidyverse)
library(Seurat)
library(scater)

# Create Shiny UI
ui <- dashboardPage(
  
  # Header
  dashboardHeader(
    title = "scRNA-Seq AI Analysis",
    tags$li(
      class = "dropdown",
      tags$a(href = "https://github.com/yourusername/scRNA-AI-App",
             "GitHub", target = "_blank")
    )
  ),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Data Loading", tabName = "data_loading", icon = icon("upload")),
      menuItem("QC & Preprocessing", tabName = "qc", icon = icon("filter")),
      menuItem("Dimensionality Reduction", tabName = "dimred", icon = icon("sitemap")),
      menuItem("Clustering", tabName = "clustering", icon = icon("object-group")),
      menuItem("Cell Type Classification", tabName = "celltype", icon = icon("dna")),
      menuItem("DE Analysis", tabName = "de_analysis", icon = icon("chart-line")),
      menuItem("AI Chatbot", tabName = "chatbot", icon = icon("robot")),
      menuItem("ML Models", tabName = "ml_models", icon = icon("brain")),
      menuItem("Reports", tabName = "reports", icon = icon("file"))
    ),
    hr(),
    p("Version 0.1.0", style = "text-align: center; font-size: 0.8em;")
  ),
  
  # Main Body
  dashboardBody(
    tabItems(
      # Dashboard Tab
      tabItem(tabName = "dashboard",
        fluidRow(
          box(title = "Welcome", width = 12,
            h3("Single Cell RNA-Seq Analysis with AI"),
            p("This application provides:"),
            tags$ul(
              tags$li("Interactive visualization of scRNA-seq data"),
              tags$li("Automated cell type classification"),
              tags$li("AI-powered chatbot for biological insights"),
              tags$li("Machine learning models for prediction"),
              tags$li("Differential expression analysis")
            ),
            br(),
            actionButton("get_started", "Get Started", class = "btn-lg btn-primary")
          )
        )
      ),
      
      # Data Loading Tab
      tabItem(tabName = "data_loading",
        h2("Data Loading"),
        fluidRow(
          box(title = "Load Data", width = 6,
            radioButtons("dataset_choice", "Select Dataset:",
              c("GSE123813 (BCC)" = "GSE123813",
                "GSE243013 (NSCLC)" = "GSE243013",
                "Upload Custom" = "custom")),
            actionButton("load_data", "Load Data", class = "btn-success"),
            br(), br(),
            uiOutput("load_status")
          ),
          box(title = "Data Summary", width = 6,
            tableOutput("data_summary_table")
          )
        )
      ),
      
      # QC Tab
      tabItem(tabName = "qc",
        h2("Quality Control"),
        fluidRow(
          box(title = "QC Metrics", width = 6,
            plotlyOutput("qc_plot")
          ),
          box(title = "Filter Settings", width = 6,
            sliderInput("min_features", "Min Features:", 
                       min = 100, max = 5000, value = 200),
            sliderInput("max_features", "Max Features:", 
                       min = 5000, max = 20000, value = 10000),
            sliderInput("max_mt", "Max Mitochondrial %:", 
                       min = 5, max = 50, value = 20),
            actionButton("apply_qc", "Apply Filters", class = "btn-info")
          )
        )
      ),
      
      # Dimensionality Reduction Tab
      tabItem(tabName = "dimred",
        h2("Dimensionality Reduction"),
        fluidRow(
          box(title = "UMAP", width = 6,
            plotlyOutput("umap_plot")
          ),
          box(title = "t-SNE", width = 6,
            plotlyOutput("tsne_plot")
          )
        ),
        fluidRow(
          box(title = "PCA Variance", width = 12,
            plotlyOutput("pca_variance_plot")
          )
        )
      ),
      
      # Clustering Tab
      tabItem(tabName = "clustering",
        h2("Clustering Analysis"),
        fluidRow(
          box(title = "Cluster Settings", width = 4,
            sliderInput("resolution", "Resolution:", 
                       min = 0.1, max = 2.0, value = 0.6),
            selectInput("clustering_algo", "Algorithm:",
                       c("Leiden" = "4", "Louvain" = "3")),
            actionButton("run_clustering", "Run Clustering", class = "btn-warning")
          ),
          box(title = "Cluster UMAP", width = 8,
            plotlyOutput("cluster_umap")
          )
        ),
        fluidRow(
          box(title = "Cluster Sizes", width = 12,
            plotlyOutput("cluster_sizes_plot")
          )
        )
      ),
      
      # Cell Type Classification Tab
      tabItem(tabName = "celltype",
        h2("Cell Type Classification"),
        fluidRow(
          box(title = "Cell Type Prediction", width = 6,
            actionButton("predict_celltypes", "Predict Cell Types", class = "btn-success"),
            br(), br(),
            uiOutput("celltype_legend")
          ),
          box(title = "Predictions", width = 6,
            plotlyOutput("celltype_umap")
          )
        ),
        fluidRow(
          box(title = "Top Marker Genes", width = 12,
            DT::dataTableOutput("marker_genes_table")
          )
        )
      ),
      
      # DE Analysis Tab
      tabItem(tabName = "de_analysis",
        h2("Differential Expression"),
        fluidRow(
          box(title = "Compare Clusters", width = 4,
            selectInput("de_cluster1", "Cluster 1:", c()),
            selectInput("de_cluster2", "Cluster 2:", c()),
            actionButton("run_de", "Run DE", class = "btn-info")
          ),
          box(title = "Volcano Plot", width = 8,
            plotlyOutput("volcano_plot")
          )
        ),
        fluidRow(
          box(title = "DE Results", width = 12,
            DT::dataTableOutput("de_results_table")
          )
        )
      ),
      
      # AI Chatbot Tab
      tabItem(tabName = "chatbot",
        h2("AI Chatbot"),
        fluidRow(
          box(title = "Chat", width = 12, height = "600px",
            div(id = "chat_history", 
              style = "height: 450px; overflow-y: auto; border: 1px solid #ccc; padding: 10px;"),
            hr(),
            fluidRow(
              column(10, textInput("chat_input", "", placeholder = "Ask about your data...")),
              column(2, actionButton("send_msg", "Send", class = "btn-primary"))
            )
          )
        ),
        fluidRow(
          box(title = "LLM Settings", width = 12,
            selectInput("llm_provider", "LLM Provider:", 
                       c("Azure OpenAI", "HuggingFace", "Local")),
            textInput("llm_api_key", "API Key:", type = "password"),
            checkboxInput("rag_enabled", "Enable RAG Context", value = TRUE)
          )
        )
      ),
      
      # ML Models Tab
      tabItem(tabName = "ml_models",
        h2("Machine Learning Models"),
        fluidRow(
          box(title = "Train Classifier", width = 6,
            actionButton("train_classifier", "Train Cell Type Classifier", 
                        class = "btn-success"),
            br(), br(),
            uiOutput("model_performance")
          ),
          box(title = "Feature Importance", width = 6,
            plotlyOutput("feature_importance_plot")
          )
        ),
        fluidRow(
          box(title = "SHAP Explanations", width = 12,
            plotlyOutput("shap_plot")
          )
        )
      ),
      
      # Reports Tab
      tabItem(tabName = "reports",
        h2("Generate Reports"),
        fluidRow(
          box(title = "Report Options", width = 12,
            checkboxGroupInput("report_sections", "Include Sections:",
              c("Summary Statistics" = "summary",
                "QC Plots" = "qc",
                "Clustering Analysis" = "clustering",
                "Cell Type Predictions" = "celltype",
                "DE Analysis" = "de",
                "Methods" = "methods")),
            actionButton("generate_report", "Generate HTML Report", 
                        class = "btn-lg btn-info")
          )
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive values to store processed data
  data_store <- reactiveValues(
    raw_seurat = NULL,
    processed_seurat = NULL,
    cell_types = NULL,
    model = NULL
  )
  
  # Data loading
  observeEvent(input$load_data, {
    shinyalert(
      title = "Loading Data",
      text = "Loading scRNA-seq data... This may take a few minutes.",
      type = "info"
    )
    
    # Placeholder - would load actual data
    output$load_status <- renderUI({
      p("✓ Data loaded successfully!", style = "color: green; font-weight: bold;")
    })
    
    output$data_summary_table <- renderTable({
      data.frame(
        Metric = c("Total Cells", "Total Genes", "Cell Types"),
        Value = c("50,000", "15,000", "12")
      )
    })
  })
  
  # Placeholder outputs
  output$qc_plot <- renderPlotly({
    plot_ly() %>%
      add_trace(x = 1:100, y = rnorm(100), type = "scatter", mode = "markers")
  })
  
  output$umap_plot <- renderPlotly({
    plot_ly() %>%
      add_trace(x = rnorm(100), y = rnorm(100), type = "scatter", mode = "markers")
  })
  
  output$tsne_plot <- renderPlotly({
    plot_ly() %>%
      add_trace(x = rnorm(100), y = rnorm(100), type = "scatter", mode = "markers")
  })
  
  # Chatbot functionality
  observeEvent(input$send_msg, {
    if (nchar(input$chat_input) > 0) {
      # Add to chat history
      message("Chat message:", input$chat_input)
      
      # Update chat display
      output$chat_history <- renderUI({
        HTML(sprintf("<p><strong>You:</strong> %s</p>", input$chat_input))
      })
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
