# Single Cell RNA-Seq Analysis with AI Chatbot & ML

A professional-grade Shiny/Golem application for interactive single-cell RNA-seq analysis, integrated with AI/LLM chatbot and machine learning models.

## 🎯 Project Overview

This project demonstrates **in-demand skills for German job market** in bioinformatics and data science:

### Core Technologies
- **R Bioinformatics**: Seurat, SingleCellExperiment, scran (scRNA-seq analysis)
- **Interactive Web**: Shiny/Golem (modular apps, reactive programming)
- **Machine Learning**: Random Forest, XGBoost, Logistic Regression, Deep Learning
- **AI/GenAI**: LLM integration (Azure OpenAI/HuggingFace), RAG, embeddings
- **Visualization**: Plotly, ggplot2 (interactive + static)
- **Data Engineering**: Large-scale processing, reproducible pipelines

### Key Skills Demonstrated
✅ Bioinformatics & Data Science  
✅ Full-Stack Web Development (R Shiny)  
✅ Machine Learning & MLOps  
✅ AI/LLM Integration  
✅ Data Visualization & Communication  
✅ Version Control & Reproducibility  

---

## 📁 Project Structure

```
scRNA_AI_App/
├── R/                          # R modules
│   ├── data_processing.R       # Load, QC, preprocessing
│   ├── ml_models.R             # Classification, clustering
│   ├── chatbot_utils.R         # LLM integration, RAG
│   └── visualization.R         # Plotting utilities
├── inst/app/
│   ├── app.R                   # Main Shiny app
│   ├── global.R                # Global config
│   └── www/                    # Static assets
├── data/
│   ├── raw/                    # Raw scRNA-seq files (GSE123813, GSE243013)
│   └── processed/              # QC'd, normalized data
├── ml_models/                  # Trained models
├── vignettes/                  # Analysis tutorials
├── tests/                      # Unit tests
├── DESCRIPTION                 # Package metadata
├── README.md                   # This file
└── .github/
    └── workflows/
        └── deploy.yml          # CI/CD pipeline
```

---

## 🚀 Quick Start

### Prerequisites
- R ≥ 4.1
- RStudio (recommended)
- Git
- Python 3.8+ (for deep learning optional)

### Installation

1. **Clone repository**
```bash
git clone https://github.com/yourusername/scRNA-AI-App.git
cd scRNA_AI_App
```

2. **Install dependencies**
```r
# Install CRAN packages
install.packages(c("shiny", "golem", "Seurat", "tidyverse", "plotly", "DT"))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scater", "scran"))

# Install from GitHub
devtools::install_github("satijalab/seurat-data")
```

3. **Configure API credentials (optional)**
```bash
# Create .env file
echo "AZURE_ENDPOINT=your_endpoint" > .env
echo "AZURE_API_KEY=your_key" >> .env
echo "HUGGINGFACE_API_KEY=your_key" >> .env
```

4. **Run the app**
```r
setwd("inst/app")
shiny::runApp()
```

---

## 📊 Features

### 1. Data Analysis
- **Load**: GSE datasets (GEO-hosted scRNA-seq)
- **QC**: Filtering, mitochondrial gene filtering, feature selection
- **Normalization**: SCTransform, log normalization
- **Dimensionality Reduction**: PCA, UMAP, t-SNE

### 2. Clustering & Cell Typing
- **Unsupervised Clustering**: Leiden/Louvain algorithms
- **Automatic Cell Type Classification**: 
  - Marker gene-based annotation
  - LLM-assisted predictions
  - Neural network classifiers
- **Differential Expression**: FindMarkers, volcano plots

### 3. AI-Powered Chatbot
- **LLM Integration**: Azure OpenAI / HuggingFace Inference
- **RAG (Retrieval-Augmented Generation)**:
  - Biological knowledge base
  - Gene expression context
  - Literature integration
- **Context-Aware Responses**: Understands your data

### 4. Machine Learning
- **Cell Type Prediction**: RF, XGBoost, Logistic Regression
- **Feature Importance**: SHAP values, permutation importance
- **Model Interpretability**: Explainable predictions
- **Deep Learning**: scVI, CellTypist (Python integration)

### 5. Advanced Analytics
- **Pathway Enrichment**: GO, KEGG pathways
- **Gene-Gene Networks**: Cell-cell communication
- **Pseudotime Analysis**: Trajectory inference
- **Batch Correction**: Harmony integration

---

## 🤖 AI Chatbot Usage Examples

### Setup
```r
# Configure LLM
llm_config <- initialize_llm(
  provider = "azure",
  endpoint = "https://yourorg.openai.azure.com/",
  api_key = Sys.getenv("AZURE_API_KEY")
)

# Enable RAG
rag_context <- create_rag_context(seurat_obj, query_genes = c("CD4", "CD8", "CD14"))
```

### Query Examples
- "What cell types are present in this dataset?"
- "Which genes differentiate cluster 0 from cluster 1?"
- "Are these T cells? What's the biological interpretation?"
- "Suggest pathway analysis for these markers"

---

## 📈 Machine Learning Pipeline

### Classification Workflow
```r
# Load and preprocess
seurat_obj <- load_and_qc_data()

# Train classifier
model <- train_celltype_classifier(
  seurat_obj,
  test_split = 0.2
)

# Get predictions
predictions <- predict_celltypes(
  new_data,
  model,
  method = "rf"
)

# Explain predictions
importance <- get_feature_importance(model, top_n = 20)
shap_vals <- generate_shap_values(model, X_test)
```

---

## 🎓 Skills for Job Search (Germany)

### Most In-Demand
1. **Bioinformatics + ML**: Seurat, differential expression, classification
2. **Full-Stack Development**: Shiny apps (frontend + backend)
3. **AI/LLM Integration**: Prompt engineering, RAG, fine-tuning
4. **Data Engineering**: ETL, processing pipelines
5. **MLOps**: Reproducibility, containerization, CI/CD

### German Job Market Focus
- Large pharma/biotech: Roche, Merck, GSK (Munich, Basel)
- Research: MPI, Helmholtz Centers
- Startups: Berlin biotech scene
- Keywords: bioinformatik, datenanalyse, KI, machine learning

### Resume Highlights
```
Skills:
• Bioinformatics: Seurat, scRNA-seq, differential expression
• Data Science: ML models, prediction, feature importance
• AI/ML Ops: LLM integration, model evaluation, SHAP
• Web Development: Shiny/Golem, reactive programming
• Data Engineering: Large-scale processing, reproducibility
• Version Control: Git, GitHub, CI/CD pipelines
```

---

## 📚 Learning Resources

### Bioinformatics
- [Seurat Tutorials](https://satijalab.org/seurat/)
- [scRNA-seq Analysis Best Practices](https://osca.bioconductor.org/)
- [Bioconductor Workflows](https://bioconductor.org/help/workflows/)

### Shiny/Golem
- [Golem Documentation](https://thinkr-open.github.io/golem/)
- [Mastering Shiny](https://mastering-shiny.org/)

### AI/LLM
- [Azure OpenAI API](https://learn.microsoft.com/en-us/azure/ai-services/openai/)
- [HuggingFace Inference API](https://huggingface.co/inference-api)
- [RAG Best Practices](https://docs.llamaindex.ai/)

### ML Interpretability
- [SHAP Documentation](https://shap.readthedocs.io/)
- [Interpretable ML Book](https://christophm.github.io/interpretable-ml-book/)

---

## 🔬 Datasets

### Included Datasets
1. **GSE123813**: Basal Cell Carcinoma (BCC) scRNA-seq
   - ~50K cells, 12 clusters
   - T cells with TCR annotation
   
2. **GSE243013**: NSCLC Immune Microenvironment
   - ~200K cells
   - Detailed immune cell annotations
   - TCR/BCR data

Both datasets are free from [GEO](https://www.ncbi.nlm.nih.gov/geo/)

---

## 🐳 Docker Deployment

```dockerfile
FROM rocker/shiny:latest

WORKDIR /home/shiny

# Install dependencies
RUN R -e "install.packages(c('shiny', 'Seurat', 'tidyverse'))"

# Copy app
COPY . .

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('inst/app')"]
```

**Build & Run**
```bash
docker build -t scrna-app .
docker run -p 3838:3838 scrna-app
```

---

## 📊 CI/CD Pipeline

GitHub Actions workflow for:
- Automated testing
- Code quality checks
- Deploy to Shiny Server
- Docker image building

See `.github/workflows/deploy.yml`

---

## 🤝 Contributing

1. Fork repository
2. Create feature branch: `git checkout -b feature/amazing-feature`
3. Commit changes: `git commit -am 'Add feature'`
4. Push to branch: `git push origin feature/amazing-feature`
5. Open Pull Request

---

## 📝 License

MIT License - see LICENSE file

---

## 🔗 Resources for Career Development

### German Biotech Companies
- **Munich**: Roche Diagnostics, Gilead, Medigene
- **Berlin**: SciLifeLab, BiotechGerm, multiple startups
- **Heidelberg**: EMBL, Germany's biotech hub
- **Frankfurt/Mainz**: Sanofi, Celgene

### Job Platforms
- LinkedIn (Tech roles in Germany)
- StepStone (German job portal)
- BioBZ (German biotech jobs)
- GitHub Jobs

### Interview Preparation
- Data science case studies
- Shiny app demo
- ML model explanation (SHAP)
- Dataset interpretation

---

## 📧 Contact & Support

- Issues: GitHub Issues
- Questions: Discussions
- Email: your.email@example.com

---

## 🎉 Acknowledgments

- Seurat team (Satija Lab, NYGENOME)
- Bioconductor community
- Single Cell Genomics Society

**Last Updated**: November 2025
**Version**: 0.1.0

---

## 🚧 Roadmap

- [ ] Deep learning integration (scVI, CellTypist)
- [ ] Spatial transcriptomics support
- [ ] Advanced RAG with vector databases
- [ ] Multi-modal analysis
- [ ] Cloud deployment templates (AWS/Azure)
- [ ] Publishable figures export
- [ ] Real-time collaboration features
