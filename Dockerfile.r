# Dockerfile for R/Shiny application
FROM rocker/shiny:4.5

LABEL maintainer="your-email@example.com"
LABEL description="Single Cell RNA-Seq Analysis with AI"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shiny', 'golem', 'DT', 'plotly', 'shinydashboard', 'shinyalert'), repos='https://cloud.r-project.org')"

RUN R -e "if (!require('BiocManager')) install.packages('BiocManager'); BiocManager::install(c('Seurat', 'SingleCellExperiment', 'scater', 'scran', 'Matrix'))"

RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'igraph', 'caret', 'randomForest', 'xgboost', 'glmnet', 'httr2', 'jsonlite', 'reticulate'), repos='https://cloud.r-project.org')"

# Set working directory
WORKDIR /home/shiny/scRNA_AI_App

# Copy application
COPY . .

# Set permissions
RUN chmod -R 755 /home/shiny

# Expose port
EXPOSE 3838

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:3838/ || exit 1

# Run Shiny app
CMD ["R", "-e", "setwd('inst/app'); shiny::runApp(port=3838, host='0.0.0.0')"]
