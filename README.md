# # Benchmarking Soft-Thresholded Compressed Sensing (ST-CS) Against Traditional Methods 
## 📂 Project Structure
```text
.
├── Data/                 # Raw proteomics data (requires download)
├── Analysis/
│ ├── Distribution.R      # Coefficient distribution visualization
│ ├── Pathway_analysis.R  # Functional enrichment analysis
│ ├── PDC000356.R         # Intrahepatic cholangiocarcinoma analysis
│ ├── PDC000446.R         # Glioblastoma analysis
│ ├── Runtime_memory.R    # Computational resource profiling
│ └── Simulation.R        # Main comparative simulation
├── Outputs/
│ ├── Figures/            # Publication-ready TIFF figures
│ └── Results/            # Analysis result tables
└── README.md
```

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
```

## 🧬 Data Preparation
- Real Data (CPTAC)
  * Download from CPTAC Data Portal (https://proteomics.cancer.gov/programs/cptac):
    + Intrahepatic cholangiocarcinoma (PDC000356):
      + iCC_NCC_Proteome.tmt10.tsv;
      + PDC_study_biospecimen_03212025_144732.csv
    + Glioblastoma (PDC000446):
      + CPTAC3_Glioblastoma_Multiforme_Confirmatory_Proteome.tmt11.tsv;
      + PDC_study_biospecimen_03232025_213615.csv
  * Store in Data/

## ▶️ Execution Workflow
### Standard Analysis
```r
# Generate coefficient distribution plot
Rscript Analysis/Distribution.R
```
```r
# Large-scale simulation study (~10 hours on 16-core system)
Rscript Analysis/Simulation.R
```
```r
# Run real data analysis (requires CPTAC data)
Rscript Analysis/PDC000356.R
Rscript Analysis/PDC000446.R
```
```r
# Perform pathway enrichment analysis
Rscript Analysis/Pathway_analysis.R
```
```r
# Benchmarking Runtime and memory usage
Rscript Analysis/Runtime_memory.R
```

## 📊 Output Specifications
| File Type	| Location	| Content | Description |
| --- | --- | --- | --- |
| TIFF Figures | Outputs/Figures/*.tiff |	Publication-quality | visualizations |
| CSV Results	| Outputs/Results/*.csv	| Quantitative analysis | metrics |


## ❓ Frequently Asked Questions
- Q: Rdonlp2 fails to install on Windows
 - A: Ensure Rtools is installed and added to PATH during installation

- Q: Implement simulation experiments with more parameters in Simulation.R
```r
# Modify in Simulation.R:
p <- 500             # Total number of features
corr <- 0.8          # Base correlation coefficient
SNR <- 10            # Signal-to-noise ratio
p_real <- 5          # Number of true signal features
n_sim <- 1000        # Number of simulations per configuration
block_size <- 50     # Size of correlation blocks
```

## 📜 Citation
```bibtex
@software{wenfudong_2025_ST-CS,
  title = {Benchmarking Soft-Thresholded Compressed Sensing (ST-CS) Against Traditional Methods},
  author = {Fudong Wen},
  year = {2025},
  url = {[https://github.com/wenfudong/ST-CS.git]},
  version = {1.0.0}
}
```

## 📄 License
MIT Licensed | Copyright © 2025 [Harbin Medical University]

