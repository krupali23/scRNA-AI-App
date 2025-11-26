# ⚡ QUICK REFERENCE CARD

## 🚀 First Time Setup (5 Minutes)

```r
# 1. Set working directory
setwd("C:/Users/krupa/Desktop/Shiny_project")

# 2. Test R version
R.version.string  # Should be 4.5.1+

# 3. Run quick start
source("QUICKSTART.R")

# 4. Launch Shiny app
shiny::runApp("inst/app")
```

## 📋 File Quick Reference

| File | Purpose | Run When |
|------|---------|----------|
| `QUICKSTART.R` | Complete example with comments | First time |
| `app.R` | Main Shiny dashboard | Launch app |
| `README.md` | Full documentation | Need help |
| `PROJECT_SUMMARY.md` | What you built & next steps | Planning |
| `ADVANCED_AI_GUIDE.md` | Cutting-edge techniques | Want to extend |
| `GITHUB_SETUP.md` | Publish to GitHub | Ready to share |
| `RESOURCES.md` | All learning resources | Need links |

## 🔧 Common Commands

### Load Data
```r
source("data_processing.R")
data <- load_gse_data("data/raw", "GSE243013")
seurat_obj <- create_seurat_object(data$counts, data$metadata)
```

### Train ML Model
```r
source("ml_models.R")
model <- train_celltype_classifier(seurat_obj, test_split = 0.2)
predictions <- predict_celltypes(new_data, model, method = "rf")
```

### Use Chatbot
```r
source("chatbot_utils.R")
llm_config <- initialize_llm(provider = "azure")
response <- query_llm("What are these cells?", "", llm_config)
```

## 🐳 Docker Commands

```bash
# Build images
docker-compose build

# Start services
docker-compose up -d

# View logs
docker-compose logs -f shiny-app

# Stop services
docker-compose down
```

## 📊 Analyze Results

### View Cluster Summary
```r
summary_table <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(n_cells = n())
```

### Find Markers
```r
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 5)
```

### Visualize
```r
# UMAP with clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Gene expression
FeaturePlot(seurat_obj, features = c("CD4", "CD8", "CD14"))

# Feature statistics
VlnPlot(seurat_obj, features = "nFeature_RNA")
```

## 🔐 Environment Setup (Optional)

```bash
# Create .env file
echo AZURE_ENDPOINT=your_endpoint > .env
echo AZURE_API_KEY=your_key >> .env
echo HUGGINGFACE_API_KEY=your_key >> .env

# Load in R
library(dotenv)
load_dot_env()
```

## 💻 GitHub Push Workflow

```bash
# Initialize
git init
git add .
git commit -m "Initial commit"

# Create remote and push
git remote add origin https://github.com/yourusername/scRNA-AI-App
git branch -M main
git push -u origin main

# Future updates
git add .
git commit -m "Add feature description"
git push
```

## 🎯 For Job Interviews

**Be ready to explain:**
1. "Walk me through your scRNA-seq pipeline"
   - Answer: Load → QC → Normalize → Reduce → Cluster → Cell type → DE

2. "How do you validate cell type predictions?"
   - Answer: Training/test split, accuracy metrics, SHAP explanations

3. "Why use multiple ML algorithms?"
   - Answer: RF for interpretability, XGBoost for accuracy, Logistic for baseline

4. "How does the chatbot work?"
   - Answer: RAG retrieves context, LLM generates response, maintains memory

5. "How would you handle 1M cells?"
   - Answer: Streaming, batch processing, dimensionality reduction

## 🌟 Impress Recruiters

1. **Show the Shiny app** (live demo is powerful!)
2. **Explain SHAP plots** (ML interpretability)
3. **Discuss data quality** (QC metrics matter)
4. **Talk about scalability** (Docker, cloud deployment)
5. **Mention AI integration** (cutting edge!)

## ⚠️ Common Issues & Fixes

| Problem | Solution |
|---------|----------|
| "Package not found" | Install: `install.packages("pkg_name")` |
| "Data file not found" | Check working directory: `getwd()` |
| "Shiny crashes" | Check R console for error messages |
| "LLM not responding" | Verify API key in .env file |
| "Out of memory" | Process in batches, reduce features |

## 📞 Help Resources

**Stuck on R?** → Check QUICKSTART.R comments  
**Shiny questions?** → See app.R tab structure  
**LLM issues?** → Review chatbot_utils.R  
**Job search?** → Read GITHUB_SETUP.md  
**Want to extend?** → See ADVANCED_AI_GUIDE.md  

## 🎓 Skills You're Demonstrating

✅ Bioinformatics (scRNA-seq, Seurat, DE)  
✅ Data Science (ML, classification)  
✅ AI/GenAI (LLM, RAG, embeddings)  
✅ Web Development (Shiny, interactive)  
✅ DevOps (Docker, Git, CI/CD)  
✅ Software Engineering (modular, documented)  

## 🚀 30-Day Plan

**Week 1**: Run QUICKSTART.R, explore Shiny app  
**Week 2**: Create GitHub repo, push code  
**Week 3**: Implement RAG chatbot  
**Week 4**: Write blog post, start job applications  

---

**Quick Win**: Run this to see everything work:
```r
source("QUICKSTART.R")
```

**Good luck! 🎉**
