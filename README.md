# # Sparse Feature Selection Framework with Soft-Thresholded Compressed Sensing (ST-CS)
## 📂 Project Structure
.
├── Data/ # Raw proteomics data (requires download)
├── Analysis/
│ ├── Distribution.R # Coefficient distribution visualization
│ ├── Pathway_analysis.R # Functional enrichment analysis
│ ├── PDC000356.R # Cholangiocarcinoma analysis
│ ├── PDC000446.R # Glioblastoma analysis
│ ├── Runtime_memory.R # Computational resource profiling
│ └── Simulation.R # Main comparative simulation
├── Outputs/
│ ├── Figures/ # Publication-ready TIFF figures
│ └── Results/ # Analysis result tables
└── README.md


## 🛠️ Environment Setup

### System Requirements
- **R ≥ 4.2.0** with LAPACK/BLAS optimization
- **Rtools** (Windows) / **Xcode** (macOS) for compiling Rdonlp2
- 8GB RAM minimum (16GB recommended for full simulations)

### Dependency Installation
```r
# Install from CRAN
install.packages(c("MASS", "Rdonlp2", "cluster", "ggplot2", "cowplot",
                   "glmnet", "spls", "pROC", "doParallel", "pryr",
                   "readr", "dplyr", "tidyr", "stringr", "tibble",
                   "matrixStats", "caret", "foreach"))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
```r

## 🧬 Data Preparation
Real Data (CPTAC)
Download from CPTAC Data Portal:
Cholangiocarcinoma: iCC_NCC_Proteome.tmt10.tsv;
                    PDC_study_biospecimen_03212025_144732.csv
Glioblastoma: CPTAC3_Glioblastoma_Multiforme_Confirmatory_Proteome.tmt11.tsv;
              PDC_study_biospecimen_03232025_213615.csv
Store in Data/
