###############################################################################
# Statistical Method Simulation Experiment
# 
# This script performs a comparative simulation study of various feature selection
# and classification methods under different experimental settings.
#
# Key Components:
# - Parallelized simulation workflow
# - Block correlation structure data generation
# - Performance evaluation metrics
# - Comparative analysis of 4 methods (STCS, HTCS, LASSO, SPLSDA)
###############################################################################

# ----------------------
# Dependencies
# ----------------------
library(doParallel)   # Parallel processing
library(foreach)      # Iteration framework
library(MASS)         # Multivariate normal distribution
library(Rdonlp2)      # Nonlinear optimization (Note: Requires special installation)
library(glmnet)       # Lasso regression
library(spls)         # Sparse partial least squares
library(cluster)      # Clustering algorithms
library(pROC)         # ROC curve analysis
library(caret)        # Performance indicator calculation
library(tidyverse)    # data processing
library(ggplot2)      # Draw graphics

# ----------------------
# Simulation Parameters
# ----------------------
p <- 500             # Total number of features
corr <- 0.8          # Base correlation coefficient
SNR <- 100           # Signal-to-noise ratio
p_real <- 5          # Number of true signal features
n_sim <- 1000        # Number of simulations per configuration
block_size <- 50     # Size of correlation blocks

# ----------------------
# Core Functions
# ----------------------

#' Generate block-structured correlation matrix
#' 
#' @param p Total number of features
#' @param block_size Size of each correlation block
#' @param corr Base correlation coefficient between features
#' @return Block-structured covariance matrix
create_block_matrix <- function(p, block_size, corr) {
  stopifnot(p %% block_size == 0)
  n_blocks <- p / block_size
  Sigma <- matrix(0, p, p)
  
  for (i in 1:n_blocks) {
    block_indices <- ((i - 1) * block_size + 1):(i * block_size)
    block_corr <- outer(1:block_size, 1:block_size, 
                        FUN = function(i, j) corr^abs(i - j))
    diag(block_corr) <- 1  # Ensure unit variance
    Sigma[block_indices, block_indices] <- block_corr
  }
  return(Sigma)
}

#' Calculate performance metrics
#' 
#' @param selected Indices of selected features
#' @param true_idx Indices of true signal features
#' @param p Total number of features
#' @param p_real Number of true signal features
#' @param pred Predicted class labels
#' @param y_true True class labels
#' @return List of performance metrics
calc_performance <- function(selected, true_idx, p, p_real, pred, y_true) {
  # Feature selection metrics
  TP <- sum(selected %in% true_idx)
  FP <- length(selected) - TP
  TN <- (p - p_real) - FP
  FN <- p_real - TP
  
  se <- TP / p_real                    # Sensitivity/Recall
  sp <- TN / (p - p_real)             # Specificity
  precision <- TP / (TP + FP)         # Precision
  fdr <- FP / (TP + FP)              # False Discovery Rate
  f1_feature <- 2 * (precision * se) / (precision + se)  # F1 Score
  
  # Matthews Correlation Coefficient
  mcc_numerator <- TP * TN - FP * FN
  mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc_feature <- ifelse(mcc_denominator == 0, 0, mcc_numerator / mcc_denominator)
  
  # Classification performance
  roc_obj <- roc(y_true, as.numeric(pred))
  auc_value <- auc(roc_obj)
  
  list(
    se = se,
    sp = sp,
    precision = precision,
    fdr = fdr,
    f1_feature = f1_feature,
    mcc_feature = mcc_feature,
    auc = auc_value
  )
}

# ----------------------
# Simulation Setup
# ----------------------
mean_x <- rep(0, p)
Sigma <- create_block_matrix(p, block_size, corr)
true_idx <- 1:p_real
beta <- rep(1, p_real)

# Initialize parallel backend
cl <- makeCluster(15)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

# ----------------------
# Main Simulation Loop
# ----------------------
results_summary <- data.frame()
start_time <- Sys.time()

for (n in seq(40, 200, 20)) {
  sim_results <- foreach(
    ii = 1:n_sim,
    .combine = rbind,
    .packages = c("MASS", "Rdonlp2", "glmnet", "cluster", "spls", "pROC", "caret"),
    .errorhandling = "remove"
  ) %dopar% {
    # Data Generation -----
    set.seed(882200 + ii)
    X_train <- mvrnorm(n, mean_x, Sigma)
    x_train <- X_train + matrix(rnorm(n * p, 0, sqrt(1 / SNR)), n, p)
    probs_train <- plogis(X_train[, true_idx] %*% beta)
    y_train <- rbinom(n, 1, probs_train)
    
    set.seed(662200 + ii)
    X_test <- mvrnorm(n, mean_x, Sigma)
    x_test <- X_test + matrix(rnorm(n * p, 0, sqrt(1 / SNR)), n, p)
    probs_test <- plogis(X_test[, true_idx] %*% beta)
    y_test <- rbinom(n, 1, probs_test)
    
    # Method Implementations -----
    
    ## 1. CS Model -----
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
    
    ## 2. STCS Method -----
    initial_medoids <- c(which.max(abs(forecast_w)), which.min(abs(forecast_w)))
    pam_clusters <- pam(abs(forecast_w), k = 2, medoids = initial_medoids)
    stcs_selected <- which(pam_clusters$clustering == which.max(pam_clusters$medoids))
    stcs_pred <- sign(x_test[, stcs_selected, drop = FALSE] %*% forecast_w[stcs_selected])
    stcs_perf <- calc_performance(stcs_selected, true_idx, p, p_real, stcs_pred, y_test * 2 - 1)
    
    ## 3. HTCS Method -----
    htcs_selected <- which(abs(forecast_w) >= 0.001)
    htcs_pred <- sign(x_test[, htcs_selected, drop = FALSE] %*% forecast_w[htcs_selected])
    htcs_perf <- calc_performance(htcs_selected, true_idx, p, p_real, htcs_pred, y_test * 2 - 1)
    
    ## 4. LASSO -----
    cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial", type.measure = "class")
    lasso_selected <- which(as.matrix(coef(cv_lasso, s = "lambda.min"))[-1] != 0)
    lasso_pred <- predict(cv_lasso, x_test, type = "class", s = "lambda.min")
    lasso_perf <- calc_performance(lasso_selected, true_idx, p, p_real, lasso_pred, y_test)
    
    ## 5. SPLSDA -----
    cv_splsda <- cv.splsda(x_train, factor(y_train), K = 1:5, 
                           eta = seq(0.1, 0.9, 0.1), fold = 5)
    splsda_model <- splsda(x_train, factor(y_train), 
                           eta = cv_splsda$eta.opt, K = cv_splsda$K.opt)
    splsda_selected <- splsda_model$A
    splsda_pred <- predict(splsda_model, x_test, fit.type = "class")
    splsda_perf <- calc_performance(splsda_selected, true_idx, p, p_real, splsda_pred, y_test)
    
    # Combine Results -----
    unlist(list(
      STCS = stcs_perf,
      HTCS = htcs_perf,
      LASSO = lasso_perf,
      SPLSDA = splsda_perf
    ))
  }
  
  # Aggregate results
  current_results <- c(n, colMeans(sim_results, na.rm = TRUE))
  results_summary <- rbind(results_summary, current_results)
}

# ----------------------
# Output Configuration
# ----------------------
metric_names <- c("se", "sp", "precision", "fdr", "f1", "mcc", "auc")
method_names <- c("STCS", "HTCS", "LASSO", "SPLSDA")
colnames(results_summary) <- c("n", 
                               paste(rep(method_names, each = length(metric_names)), 
                                     metric_names, 
                                     sep = "_"))

# SNR-specific output files
output_file <- file.path("Outputs", "Results", paste0("Simulation_results_SNR_", SNR, ".csv"))
write.csv(results_summary, output_file, row.names = FALSE)

# ----------------------
# Visualization 
# ----------------------

# convert it to a long format
long_data <- results_summary %>%
  pivot_longer(
    cols = -n,
    names_to = c("Method", "Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  filter(Metric %in% c("se", "sp", "precision", "fdr", "f1", "mcc")) %>%
  mutate(
    Method = recode(Method,
                    STCS = "ST-CS",
                    HTCS = "HT-CS",
                    LASSO = "LASSO",
                    SPLSDA = "SPLSDA"
    ),
    Metric = factor(Metric,
                    levels = c("se", "sp", "precision", "fdr", "f1", "mcc"),
                    labels = c("(a) Sensitivity", "(b) Specificity", 
                               "(c) Precision", "(d) FDR", "(e) F1 Score", "(f) MCC")
    )
  )

# Customize Color Palette
method_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")

# Draw graphics
ggplot(long_data, aes(x = n, y = Value, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~Metric, nrow = 2, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(x = "Sample Size", y = "Metric Value") +
  scale_color_manual(values = method_colors) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10, hjust = 0),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 1, title.position = "top"))

figure_file <- file.path("Outputs", "Figures", paste0("Comparison_SNR_", SNR, ".tiff"))
ggsave(figure_file, width = 12, height = 8, dpi = 300)

cat("Completed SNR =", SNR, "Runtime:", Sys.time() - start_time, "\n")

