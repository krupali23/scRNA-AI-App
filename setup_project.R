# Setup AI-Enhanced Single Cell Analysis Shiny App
# This script creates a professional Golem app structure

# Create project directories
if (!dir.exists("scRNA_AI_App")) {
  dir.create("scRNA_AI_App", recursive = TRUE)
}

setwd("scRNA_AI_App")

# Create essential project structure
dirs <- c(
  "R",
  "inst",
  "inst/app",
  "data",
  "data-raw",
  "man",
  "tests",
  "vignettes",
  "data/processed",
  "data/raw",
  "ml_models",
  "shiny_modules"
)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Create DESCRIPTION file
description <- '
Package: scrnaAIApp
Title: AI-Enhanced Single Cell RNA-Seq Analysis
Version: 0.1.0
Authors@R: person("Your Name", "Your Email", role = c("aut", "cre"))
Description: Interactive Shiny/Golem application for single cell RNA-seq analysis with integrated AI chatbot and ML models.
License: MIT
Imports:
    shiny,
    golem,
    DT,
    ggplot2,
    dplyr,
    Seurat,
    SingleCellExperiment,
    scater,
    scran,
    plotly,
    httr,
    jsonlite,
    reticulate,
    caret,
    igraph,
    tidyverse
Suggests:
    testthat,
    knitr,
    rmarkdown
VignetteBuilder: knitr
'

writeLines(trimws(description), "DESCRIPTION")

# Create global.R
global_content <- '
# Global configuration for scRNA-seq AI analysis app
library(shiny)
library(Seurat)
library(tidyverse)
library(plotly)
library(DT)

# Set data paths
DATA_PATH <- here::here("data")
RAW_DATA_PATH <- file.path(DATA_PATH, "raw")
PROCESSED_DATA_PATH <- file.path(DATA_PATH, "processed")
MODELS_PATH <- here::here("ml_models")

# API configuration for AI chatbot
AI_API_ENDPOINT <- Sys.getenv("AI_API_ENDPOINT", "")
AI_API_KEY <- Sys.getenv("AI_API_KEY", "")

# Create directories if they don't exist
dir.create(PROCESSED_DATA_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(MODELS_PATH, showWarnings = FALSE, recursive = TRUE)

# Source utility functions
source("R/data_processing.R")
source("R/ml_models.R")
source("R/chatbot_utils.R")
source("R/visualization.R")
'

writeLines(trimws(global_content), "inst/app/global.R")

cat("✓ Project structure created successfully!\n")
cat("✓ Location: ", getwd(), "\n")
cat("✓ Navigate to this directory and run setup_golem.R next\n")
