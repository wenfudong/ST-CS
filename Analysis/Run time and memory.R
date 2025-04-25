# ------------------------------------------------------------------------------
# Code Overview: High-Dimensional Constrained Optimization Benchmark
#
# This script evaluates computational performance (runtime & memory usage) of 
# constrained nonlinear optimization under increasing feature dimensions (p). 
#
# Key components:
# 1. Synthetic Data Generation
# 2. Optimization: Solves L1/L2-constrained maximum margin problem using Rdonlp2
# 3. Performance Profiling
# ------------------------------------------------------------------------------

# ----------------------
# Dependencies
# ----------------------
library(MASS)         # For multivariate normal distribution
library(Rdonlp2)      # For nonlinear constrained optimization
library(pryr)         # For memory usage tracking

# ----------------------
# Simulation Parameters
# ----------------------
corr <- 0.8          # Base correlation coefficient between features
SNR <- 100           # Signal-to-noise ratio for data generation
p_real <- 5          # Number of true causal features
block_size <- 50     # Size of correlation blocks in design matrix
n <- 200             # Training sample size

# ----------------------
# Core Functions
# ----------------------

#' Generate block-structured covariance matrix
#' 
#' @param p Total number of features (must be multiple of block_size)
#' @param block_size Number of features per correlation block
#' @param corr Base correlation coefficient within blocks
#' @return Block-diagonal covariance matrix with AR(1) structure within blocks
create_block_matrix <- function(p, block_size, corr) {
  stopifnot(p %% block_size == 0)
  n_blocks <- p / block_size
  Sigma <- matrix(0, p, p)
  
  # Create AR(1) correlation structure within each block
  for (i in 1:n_blocks) {
    block_indices <- ((i - 1) * block_size + 1):(i * block_size)
    block_corr <- outer(1:block_size, 1:block_size, 
                        FUN = function(i, j) corr^abs(i - j))
    diag(block_corr) <- 1  # Ensure unit variance
    Sigma[block_indices, block_indices] <- block_corr
  }
  return(Sigma)
}

# ----------------------
# Initialize Results Storage
# ----------------------
results <- data.frame(
  p = integer(),      # Feature dimension
  time_sec = numeric(),  # Computation time (seconds)
  memory_mb = numeric()  # Memory usage (MB)
)

# ----------------------
# Simulation Loop
# ----------------------
for (p in seq(2000, 12000, 2000)) {
  # ----------------------
  # Data Generation
  # ----------------------
  # Create covariance matrix with block structure
  mean_x <- rep(0, p)
  Sigma <- create_block_matrix(p, block_size, corr)
  
  # True signal indices and coefficients
  true_idx <- 1:p_real
  beta <- rep(1, p_real)
  
  # Generate training data with additive noise
  set.seed(882200 + p)  # For reproducibility
  X_train <- mvrnorm(n, mean_x, Sigma)  # Latent features
  x_train <- X_train + matrix(rnorm(n * p, 0, sqrt(1/SNR)), n, p)  # Observed features
  probs_train <- plogis(X_train[, true_idx] %*% beta)  # True probabilities
  y_train <- rbinom(n, 1, probs_train)  # Binary responses
  
  # ----------------------
  # Resource Monitoring
  # ----------------------
  gc()  # Garbage collection to ensure accurate memory measurement
  
  # Baseline memory usage
  mem_before <- mem_used() 
  
  # Start timing
  start_time <- Sys.time()
  
  # ----------------------
  # Constrained Optimization Model
  # ----------------------
  # Objective function: maximize margin
  objective_fn <- function(w) -(y_train * 2 - 1) %*% x_train %*% w
  
  # Initial parameter values
  w_init <- rep(0, p)
  
  # Constraints: L1 and L2 norm bounds
  constraints <- list(
    function(w) sum(abs(w)) - sqrt(1),  # L1 constraint
    function(w) sqrt(sum(w^2)) - 1      # L2 constraint
  )
  
  # Solve constrained optimization problem
  cs_fit <- donlp2(
    w_init, objective_fn,
    par.upper = rep(1, p), par.lower = rep(-1, p),
    nlin.upper = rep(0, 2), nlin.lower = rep(-Inf, 2),
    nlin = constraints
  )
  
  # ----------------------
  # Final Measurements
  # ----------------------
  time_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  mem_after <- mem_used()
  
  # ----------------------
  # Record Results
  # ----------------------
  results <- rbind(results, data.frame(
    p = p,
    time_sec = round(as.numeric(time_elapsed), 2),
    memory_mb = round((mem_after - mem_before)/1024^2, 2)  # Convert bytes to MB
  ))
  
  # Progress update
  message(sprintf("Completed p=%d: Time=%.1fs, Memory=%.2fMB",
                  p, results$time_sec[nrow(results)], results$memory_mb[nrow(results)]))
}

# ----------------------
# Output Results
# ----------------------
print(results)
write.csv(results, 
          "Outputs/Results/Optimization_benchmark.csv", 
          row.names = FALSE)  # Save results for analysis