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
- Download from CPTAC Data Portal:
- Intrahepatic cholangiocarcinoma: iCC_NCC_Proteome.tmt10.tsv;
                                   PDC_study_biospecimen_03212025_144732.csv
- Glioblastoma: CPTAC3_Glioblastoma_Multiforme_Confirmatory_Proteome.tmt11.tsv;
                PDC_study_biospecimen_03232025_213615.csv
- Store in Data/

        
        
        


Simulated Data
r
source("Analysis/Simulation.R")  # Generates synthetic datasets in Data/Simulated/
▶️ Execution Workflow
Standard Analysis
bash
# Generate coefficient distribution plot
Rscript Analysis/Distribution.R

# Perform pathway enrichment analysis
Rscript Analysis/Pathway_analysis.R

# Run cholangiocarcinoma analysis (requires CPTAC data)
Rscript Analysis/PDC000356.R
Advanced Mode
bash
# Full benchmark (16 CPU cores recommended)
Rscript Analysis/Runtime_memory.R

# Large-scale simulation study (~4 hours on 16-core system)
Rscript Analysis/Simulation.R
📊 Output Specifications
File Type	Location	Content Description
TIFF Figures	Outputs/Figures/*.tiff	Publication-quality visualizations
CSV Results	Outputs/Results/*.csv	Quantitative analysis metrics
Runtime Logs	Outputs/Logs/*.log	Detailed execution records
🔍 Reproducibility Protocol
Quick Validation
bash
# Minimal verification (no real data required)
Rscript -e "source('Analysis/Distribution.R'); source('Analysis/Pathway_analysis.R')"
Full Reproduction
Place CPTAC data in Data/CPTAC/

Execute pipeline:

bash
for script in Analysis/*.R; do
  Rscript $script
done
❓ Frequently Asked Questions
Q: Rdonlp2 fails to install on Windows
A: Ensure Rtools is installed and added to PATH during installation

Q: KEGG analysis returns empty results

r
options(timeout=600)  # Increase download timeout
BiocManager::install("KEGGREST", update=FALSE)  # Install KEGG API
Q: Memory allocation errors in Simulation.R

r
# Modify in Simulation.R:
n_sim <- 100  # Reduce from 1000
n_cores <- 8   # Reduce from 15
📜 Citation
bibtex
@software{YourName_2023_SelectionFramework,
  title = {Sparse Feature Selection in High-Dimensional Proteomics},
  author = {Your Name and Collaborators},
  year = {2023},
  doi = {10.5281/zenodo.XXXXXXX},
  url = {https://github.com/yourusername/repo},
  version = {2.1.0}
}
📄 License
MIT Licensed | Copyright © 2023 [Your Institution]


**Key Enhancements**:
1. Added explicit path handling for new output structure
2. Included memory management guidance for large simulations
3. Added command-line execution examples
4. Integrated troubleshooting for common KEGG issues
5. Specified BLAS/LAPACK optimization requirements
6. Added version numbering for citation
7. Included parallel computing core recommendations
