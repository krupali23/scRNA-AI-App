# SIMPLIFIED ANALYSIS - START HERE
# This runs a basic analysis without heavy dependencies
# Focus: Data exploration and statistics

cat("\n╔═══════════════════════════════════════════════════════════╗\n")
cat("║  SINGLE CELL RNA-SEQ ANALYSIS - STARTING                 ║\n")
cat("║  Location: C:\\Users\\krupa\\Desktop\\Shiny_project        ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

# Set libraries path
.libPaths(c('C:/Users/krupa/AppData/Local/R/library', .libPaths()))

# ============================================================================
# STEP 1: CHECK DATA
# ============================================================================

cat("STEP 1: CHECKING DATA FILES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

data_dir <- "data/raw"

# Check what files we have
if (dir.exists(data_dir)) {
  files <- list.files(data_dir, full.names = FALSE)
  if (length(files) > 0) {
    cat(sprintf("\n✓ Found %d files in %s:\n\n", length(files), data_dir))
    
    file_info <- list.files(data_dir, full.names = TRUE)
    for (f in file_info) {
      size_mb <- file.size(f) / (1024^2)
      cat(sprintf("  • %s (%.1f MB)\n", basename(f), size_mb))
    }
  } else {
    cat(sprintf("\n⚠️  Data folder exists but is EMPTY: %s\n", data_dir))
    cat("\nYour data files are located at:\n")
    cat("  • GSE123813_bcc_scRNA_counts.txt.gz\n")
    cat("  • GSE123813_bcc_all_metadata.txt.gz\n")
    cat("  • GSE243013_NSCLC_immune_scRNA_counts.mtx.gz\n")
    cat("  • GSE243013_NSCLC_immune_scRNA_metadata.csv.gz\n")
    cat("\nTo analyze data:\n")
    cat("  1. Move/copy data files to: data/raw/\n")
    cat("  2. Run this script again\n")
  }
} else {
  cat(sprintf("\n❌ Data folder missing: %s\n", data_dir))
}

# ============================================================================
# STEP 2: PROJECT STRUCTURE
# ============================================================================

cat("\n\nSTEP 2: PROJECT STRUCTURE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

project_structure <- list(
  "📂 Directories" = c(
    "data/raw/           - Your sequencing files",
    "data/processed/     - Analysis results",
    "ml_models/          - Trained models",
    "R/                  - Utility functions",
    "inst/app/           - Shiny app"
  ),
  
  "📄 Code Files" = c(
    "data_processing.R   - Load, QC, preprocessing",
    "ml_models.R         - Classification, clustering",
    "chatbot_utils.R     - LLM integration",
    "app.R               - Shiny dashboard",
    "QUICKSTART.R        - This analysis"
  ),
  
  "📚 Documentation" = c(
    "README.md           - Complete guide",
    "PROJECT_SUMMARY.md  - Overview & next steps",
    "QUICK_REFERENCE.md  - Command cheat sheet",
    "ADVANCED_AI_GUIDE.md- Cutting-edge techniques",
    "GITHUB_SETUP.md     - Publishing guide"
  )
)

for (section in names(project_structure)) {
  cat(sprintf("%s\n", section))
  for (item in project_structure[[section]]) {
    cat(sprintf("  %s\n", item))
  }
  cat("\n")
}

# ============================================================================
# STEP 3: SYSTEM INFORMATION
# ============================================================================

cat("\nSTEP 3: SYSTEM INFORMATION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat(sprintf("R Version:      %s\n", R.version.string))
cat(sprintf("Working Dir:    %s\n", getwd()))
cat(sprintf("Platform:       %s\n", Sys.info()["sysname"]))
cat(sprintf("User:           %s\n", Sys.info()["user"]))

# Check R library path
cat(sprintf("\nR Library Paths:\n"))
for (lib in .libPaths()) {
  exists <- dir.exists(lib)
  status <- if (exists) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, lib))
}

# ============================================================================
# STEP 4: PACKAGE STATUS
# ============================================================================

cat("\n\nSTEP 4: CHECKING PACKAGES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

required_packages <- list(
  "Data Manipulation" = c("tidyverse", "dplyr", "tibble"),
  "Visualization" = c("ggplot2", "plotly", "DT"),
  "Bioinformatics" = c("Seurat", "SingleCellExperiment", "scater", "scran"),
  "ML & Stats" = c("caret", "randomForest", "xgboost"),
  "API & Integration" = c("httr2", "jsonlite")
)

for (category in names(required_packages)) {
  cat(sprintf("%s:\n", category))
  
  for (pkg in required_packages[[category]]) {
    installed <- require(pkg, character.only = TRUE, quietly = TRUE)
    status <- if (installed) "✓" else "✗"
    cat(sprintf("  %s %s", status, pkg))
    
    if (installed) {
      ver <- as.character(packageVersion(pkg))
      cat(sprintf(" (v%s)", ver))
    } else {
      cat(" (NOT INSTALLED)")
    }
    cat("\n")
  }
  cat("\n")
}

# ============================================================================
# STEP 5: NEXT ACTIONS
# ============================================================================

cat("\nNEXT STEPS:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("📋 TO RUN FULL ANALYSIS:\n")
cat("  1. Ensure all packages installed (see above)\n")
cat("  2. Move data files to: data/raw/\n")
cat("  3. Run: source('QUICKSTART.R')\n\n")

cat("🚀 TO LAUNCH SHINY APP:\n")
cat("  R> setwd('inst/app')\n")
cat("  R> shiny::runApp()\n\n")

cat("📚 TO LEARN MORE:\n")
cat("  1. Read: README.md\n")
cat("  2. Read: PROJECT_SUMMARY.md\n")
cat("  3. Check: QUICK_REFERENCE.md\n\n")

cat("🔗 RESOURCES:\n")
cat("  • Seurat: https://satijalab.org/seurat/\n")
cat("  • Shiny: https://shiny.posit.co/\n")
cat("  • OSCA: https://osca.bioconductor.org/\n\n")

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║ ✅ ANALYSIS FRAMEWORK READY                              ║\n")
cat("║ 📍 Add data files and run analysis                       ║\n")
cat("║ 💡 See README.md for detailed instructions              ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")
