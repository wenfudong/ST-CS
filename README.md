# ST-CS
## 📂 Project Structure
.
├── Analysis/
│ ├── Distribution.R # Coefficient distribution visualization
│ ├── Pathway_analysis.R # GO/KEGG enrichment analysis
│ ├── PDC000356.R # Intrahepatic cholangiocarcinoma analysis
│ ├── PDC000446.R # Glioblastoma analysis
│ ├── Runtime_memory.R # Computational resource evaluation
│ ├── Simulation.R # Main simulation experiment
│ └── Visualization of simulations.R # Visualization of simulation experiment results
├── Outputs/ # Generated results
│ ├── Figures/ # Publication-quality figures
│ └── Results/ # Analysis result tables
└── README.md # Project documentation

## 🛠️ Environment Setup

### System Requirements
- R ≥ 4.2.0
- Rtools (Windows users)
- 8GB+ RAM (16GB recommended for full simulations)

### Dependency Installation
```r
# CRAN packages
install.packages(c("MASS", "Rdonlp2", "cluster", "ggplot2", "cowplot",
                   "glmnet", "spls", "pROC", "doParallel", "pryr"))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
