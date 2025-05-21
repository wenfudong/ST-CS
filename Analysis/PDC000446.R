# -----------------------------------------------------------------------
# Cancer Proteomics Analysis Pipeline (PDC000446)
#
# This script performs comparative analysis of feature selection methods
# for cancer classification using proteomics data. Includes four methods:
# 1. Soft Thresholding CS (STCS)
# 2. Hard Thresholding CS (HTCS)
# 3. LASSO Regression
# 4. Sparse PLS-DA
#
# Output: CSV file with comparative performance metrics
# ----------------------------------------------------------------------

# Load required packages
library(readr)       # Efficient data reading
library(dplyr)       # Data manipulation
library(tidyr)       # Data tidying
library(stringr)     # String manipulation
library(tibble)      # Modern data frames
library(matrixStats) # Matrix statistics
library(Rdonlp2)     # Nonlinear optimization (Note: Requires special installation)
library(glmnet)      # Regularized regression
library(spls)        # Sparse partial least squares
library(cluster)     # Clustering algorithms
library(pROC)        # ROC curve analysis

# ----------------------------
# Data Loading
# ----------------------------

# Read proteomic data
proteomic_data <- read_tsv("Data/CPTAC3_Glioblastoma_Multiforme_Confirmatory_Proteome.tmt11.tsv")

# Read clinical data
clinical_data <- read_csv("Data/PDC_study_biospecimen_03232025_213615.csv")

# ----------------------------
# Data Preprocessing
# ----------------------------

# Process proteomics data
log_intensity_data <- proteomic_data %>%
  column_to_rownames(var = "Gene") %>%          # Set gene names as row names
  select(contains("Log ") & !contains("Unshared")) %>%     # Select log-transformed features
  filter(!rownames(.) %in% c("Mean", "Median", "StdDev"))  # Remove summary rows

# Clean column names by extracting sample identifiers
colnames(log_intensity_data) <- str_remove_all(
  colnames(log_intensity_data), 
  " Log Ratio"
  )

# Filter clinical samples
clinical_samples <- clinical_data %>% 
  filter(`Sample Type` %in% c("Primary Tumor", "Recurrent Tumor"),
         `Aliquot Status` == "Qualified")

# Create binary outcome vector (0: Primary, 1: Recurrent)
y <- ifelse(clinical_samples$`Sample Type` == "Primary Tumor", 0, 1)

# Align proteomic data with clinical samples
sample_ids <- clinical_samples$`Aliquot Submitter ID`
aligned_data <- log_intensity_data[, sample_ids, drop = FALSE] %>% 
  na.omit() %>%     # Remove missing values
  t()               # Transpose to samples x features

# ----------------------------
# Feature Selection
# ----------------------------

# Remove low-variance features (keep top 20% highest variance)
feature_variances <- matrixStats::rowVars(t(aligned_data))
variance_threshold <- quantile(feature_variances, 0.8)
selected_features <- which(feature_variances > variance_threshold)

# Standardize selected features (mean=0, sd=1)
processed_data <- scale(aligned_data[, selected_features])

# ----------------------------
# Cross-Validation Setup
# ----------------------------

set.seed(12345)  # Reproducibility seed
n_folds <- 5     # Number of cross-validation folds
fold_ids <- sample(rep(1:n_folds, length.out = nrow(processed_data)))

# Initialize result containers
results <- list(
  stcs = list(auc = numeric(n_folds), n_features = numeric(n_folds)),
  htcs = list(auc = numeric(n_folds), n_features = numeric(n_folds)),
  lasso = list(auc = numeric(n_folds), n_features = numeric(n_folds)),
  splsda = list(auc = numeric(n_folds), n_features = numeric(n_folds))
)

# Initialize ST-CS feature collection container
all_stcs_features <- character()

# ----------------------------
# Main Cross-Validation Loop
# ----------------------------

for (fold in 1:n_folds) {
  # Data partitioning
  test_indices <- which(fold_ids == fold)
  x_train <- processed_data[-test_indices, ]
  y_train <- y[-test_indices]
  x_test <- processed_data[test_indices, ]
  y_test <- y[test_indices]
  
  # ------------------------
  # 1. Compressive Sensing (CS) Model
  # ------------------------
  # Quadratic programming implementation
  objective_fn <- function(w) -(y_train * 2 - 1) %*% x_train %*% w
  w_init <- rep(0, ncol(x_train))
  constraints <- list(
    function(w) sum(abs(w)) - sqrt(1),  # L1 constraint
    function(w) sqrt(sum(w^2)) - 1      # L2 constraint
  )
  
  cs_fit <- donlp2(
    w_init, objective_fn,
    par.upper = rep(1, ncol(x_train)), par.lower = rep(-1, ncol(x_train)),
    nlin.upper = rep(0, 2), nlin.lower = rep(-Inf, 2),
    nlin = constraints
  )
  forecast_w <- cs_fit$par
  
  # ------------------------
  # 2. STCS Method (Soft Thresholding)
  # ------------------------
  initial_medoids <- c(which.max(abs(forecast_w)), which.min(abs(forecast_w)))
  pam_clusters <- pam(abs(forecast_w), k = 2, medoids = initial_medoids)
  stcs_selected <- which(pam_clusters$clustering == which.max(pam_clusters$medoids))
  stcs_pred <- sign(x_test[, stcs_selected, drop = FALSE] %*% forecast_w[stcs_selected])
  results$stcs$auc[fold] <- pROC::roc(y_test, as.numeric(stcs_pred))$auc
  results$stcs$n_features[fold] <- length(stcs_selected)
  
  # Record the selected feature names
  stcs_feature_names <- colnames(processed_data)[stcs_selected]
  all_stcs_features <- c(all_stcs_features, stcs_feature_names)
  
  # ------------------------
  # 3. HTCS Method (Hard Thresholding)
  # ------------------------
  htcs_selected <- which(abs(forecast_w) >= 0.001)
  htcs_pred <- sign(x_test[, htcs_selected, drop = FALSE] %*% forecast_w[htcs_selected])
  results$htcs$auc[fold] <- pROC::roc(y_test, as.numeric(htcs_pred))$auc
  results$htcs$n_features[fold] <- length(htcs_selected)
  
  # ------------------------
  # 4. Regularized Regression (LASSO)
  # ------------------------
  cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial", type.measure = "class")
  lasso_selected <- which(as.matrix(coef(cv_lasso, s = "lambda.min"))[-1] != 0)
  lasso_pred <- predict(cv_lasso, x_test, type = "class", s = "lambda.min")
  results$lasso$auc[fold] <- pROC::roc(y_test, as.numeric(lasso_pred))$auc
  results$lasso$n_features[fold] <- length(lasso_selected)
  
  # ------------------------
  # 5. Sparse PLS-DA
  # ------------------------
  cv_splsda <- cv.splsda(x_train, factor(y_train), K = 1:5, 
                         eta = seq(0.1, 0.9, 0.1), fold = 5)
  splsda_model <- splsda(x_train, factor(y_train), 
                         eta = cv_splsda$eta.opt, K = cv_splsda$K.opt)
  splsda_selected <- splsda_model$A
  splsda_pred <- predict(splsda_model, x_test, fit.type = "class")
  results$splsda$auc[fold] <- pROC::roc(y_test, as.numeric(splsda_pred))$auc
  results$splsda$n_features[fold] <- length(splsda_selected)
}

# ----------------------------
# Results Export to CSV
# ----------------------------

# Create summary data frame
result_summary <- data.frame(
  Method = c("STCS", "HTCS", "LASSO", "SPLS-DA"),
  Mean_AUC = sapply(results, function(x) round(mean(x$auc), 4)),
  SD_AUC = sapply(results, function(x) round(sd(x$auc), 4)),
  Mean_Features = sapply(results, function(x) round(mean(x$n_features))),
  SD_Features = sapply(results, function(x) round(sd(x$n_features)))
)

# Write to CSV file
write_csv(result_summary, "Outputs/Results/method_comparison_PDC000446.csv")

# ----------------------------
# Feature frequency analysis module
# ----------------------------
# Create a frequency data box and sort it
feature_frequency <- data.frame(
  Protein = names(table(all_stcs_features)),
  Frequency = as.numeric(table(all_stcs_features))
) %>%
  arrange(desc(Frequency)) 

# Save feature frequency results
write_csv(feature_frequency, "Outputs/Results/STCS_feature_PDC000446.csv")

