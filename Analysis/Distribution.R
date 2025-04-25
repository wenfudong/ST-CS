# ==============================================================================
# Single Simulation of Sparse Feature Selection with Clustering Visualization
# ==============================================================================

# -------------------------------
# 1. Environment Configuration
# -------------------------------
# 1.1 Set reproducibility seed
set.seed(2023)

# 1.2 Load libraries
library(MASS)         # For mvrnorm
library(Rdonlp2)      # For nonlinear optimization (Note: Requires Rtools on Windows)
library(cluster)      # For pam clustering
library(ggplot2)      # For visualization
library(cowplot)      # For plot arrangement

# -------------------------------
# 2. Core Simulation Parameters
# -------------------------------
SIM_PARAMS <- list(
  p = 500,            # Total features
  n = 200,            # Sample size
  block_size = 50,    # Correlation block size
  corr = 0.8,         # Base correlation
  SNR = 100,          # Signal-to-noise ratio
  p_real = 5          # True biomarkers
)

# -------------------------------
# 3. Core Functions
# -------------------------------

#' Generate block-structured covariance matrix
#' @description Creates autoregressive block correlation structure
generate_block_matrix <- function(p, block_size, base_corr) {
  validate_matrix_params(p, block_size)
  
  n_blocks <- p / block_size
  Sigma <- matrix(0, p, p)
  
  for (i in 1:n_blocks) {
    indices <- ((i-1)*block_size+1):(i*block_size)
    block <- outer(1:block_size, 1:block_size, 
                   function(x,y) base_corr^abs(x-y))
    Sigma[indices, indices] <- block
  }
  diag(Sigma) <- 1 # Ensure unit variance
  Sigma
}

validate_matrix_params <- function(p, block_size) {
  if (p %% block_size != 0) stop("p must be divisible by block_size")
}

#' Generate simulated proteomics data
generate_sim_data <- function(n, p, p_real, Sigma, SNR) {
  X_clean <- mvrnorm(n, rep(0, p), Sigma)
  noise <- matrix(rnorm(n*p, 0, sqrt(1/SNR)), n, p)
  X_noisy <- X_clean + noise
  y <- rbinom(n, 1, plogis(X_clean[,1:p_real] %*% rep(1, p_real)))
  list(X = X_noisy, y = y)
}

#' Fit compressed sensing model
fit_cs_model <- function(X, y) {
  obj_fn <- function(w) -(2*y - 1) %*% X %*% w
  constraints <- list(
    function(w) sum(abs(w)) - sqrt(1),  # L1
    function(w) sqrt(sum(w^2)) - 1      # L2
  )
  
  donlp2(par = rep(0, ncol(X)),
         fn = obj_fn,
         par.upper = rep(1, ncol(X)),
         par.lower = rep(-1, ncol(X)),
         nlin.upper = rep(0, 2), 
         nlin.lower = rep(-Inf, 2),
         nlin = constraints)
}

# -------------------------------
# 4. Simulation Workflow
# -------------------------------

# 4.1 Generate covariance matrix
Sigma <- generate_block_matrix(SIM_PARAMS$p, SIM_PARAMS$block_size, SIM_PARAMS$corr)

# 4.2 Generate simulated data
sim_data <- generate_sim_data(SIM_PARAMS$n, SIM_PARAMS$p, SIM_PARAMS$p_real,
                              Sigma, SIM_PARAMS$SNR)

# 4.3 Fit compressed sensing model
cs_result <- fit_cs_model(sim_data$X, sim_data$y)

# 4.4 Extract coefficients
coefficients <- cs_result$par
abs_coeff <- abs(coefficients)

# -------------------------------
# 5. Clustering Analysis
# -------------------------------

# 5.1 Perform K-Medoids clustering
initial_medoids <- c(which.max(abs_coeff), which.min(abs_coeff))
kmed <- pam(abs_coeff, k = 2, medoids = initial_medoids)
cluster_labels <- factor(kmed$clustering, 
                         levels = c(1,2),
                         labels = c("Biomarkers", "Noise"))

# 5.2 Create analysis dataframe
analysis_df <- data.frame(
  coefficient = abs_coeff,
  cluster = cluster_labels
)

# -------------------------------
# 6. Visualization
# -------------------------------

#' Create density plot for cluster visualization
create_cluster_plot <- function(data, cluster_name, color) {
  ggplot(data, aes(x = coefficient)) +
    geom_density(adjust = 2,
                 color = color,
                 linewidth = 0.5) +
    geom_area(stat = "density",
              adjust = 2,
              fill = color,
              alpha = 0.2) +
    labs(title = cluster_name,
         x = expression("|"*omega*"*|"),
         y = "Density") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
}

# 6.2 Generate individual plots
noise_plot <- create_cluster_plot(
  subset(analysis_df, cluster == "Noise"),
  "Noise Cluster",
  "#191970"
)

biomarker_plot <- create_cluster_plot(
  subset(analysis_df, cluster == "Biomarkers"),
  "Biomarker Cluster",
  "#FF7F50"
)

# 6.3 Combine plots
final_plot <- plot_grid(
  noise_plot,
  biomarker_plot,
  ncol = 2,
  labels = "AUTO",
  align = "hv"
)

# 6.4 Save output
ggsave("Outputs/Figures/Coefficient_distribution.tiff",
       final_plot,
       device = "tiff",
       width = 10,
       height = 5,
       dpi = 300)

# 6.5 Display plot
print(final_plot)