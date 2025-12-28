##########################################
# Real Dataset 1: Riboflavin Production Dataset (top 100 genes)
# Goal: Evaluate SCAD, MCP, LOG, Hybrid, and Ensemble methods on noisy + missing data

# Load required packages
library(hdi)
library(MASS)
library(Matrix)
library(Metrics)
library(dplyr)

# Data: Top 100 genes from Riboflavin 
# Prepare Data: Top 100 genes from Riboflavin 
data(riboflavin)
X_full <- scale(riboflavin$x)
top100 <- order(apply(X_full, 2, var), decreasing = TRUE)[1:100]
X <- X_full[, top100]
n <- nrow(X)
p <- ncol(X)

#  Define parameters 
lambda <- 0.3
missing_rate <- 0.2
outlier_rate <- 0.1
contam_sd <- 25
runs <- 50 # run 50 times

#  Real data experiment runner for any penalty 
run_real_experiment <- function(penalty = "scad", runs = 10) {
  results <- vector("list", runs)
  overall_start <- Sys.time()
  for (i in 1:runs) {
    #cat("Run", i, "with", penalty, "\n")
    
    # Simulate sparse beta and y
    true_beta <- rep(0, p)
    #nonzero_idx <- sample(1:p, 10)
    #true_beta[nonzero_idx] <- runif(10, -2, 2)
    nonzero_count <- min(10, p) # fit for the p<10 case dataset
    nonzero_idx <- sample(1:p, nonzero_count)
    true_beta[nonzero_idx] <- runif(nonzero_count, -2, 2)
    y <- X %*% true_beta + rnorm(n)
    y <- scale(y)
    
    # Inject outliers
    outlier_idx <- sample(1:n, floor(outlier_rate * n))
    y[outlier_idx] <- y[outlier_idx] + rnorm(length(outlier_idx), 0, contam_sd)
    
    # Inject missing values
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * missing_rate)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    
    # Impute with mean
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    
    # EM iterations
    beta <- rep(0.01, p)
    start_time <- Sys.time()
    for (em in 1:10) {
      X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
      colnames(X_df) <- paste0("V", 1:p)
      for (j in 1:p) {
        obs <- !is.na(X_miss[, j])
        if (sum(obs) < n) {
          fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))), data = X_df[obs, ]), silent = TRUE)
          if (!inherits(fit, "try-error")) {
            X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
          }
        }
      }
      X_imp <- as.matrix(X_df)
      
      if (penalty == "hybrid") {
        beta <- coordinate_descent_hybrid(X_imp, y, beta, lambda)
      } else {
        beta <- coordinate_descent(X_imp, y, beta, lambda, penalty)
      }
    }
    runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    pred_nonzero <- which(abs(beta) > 1e-4)
    
    # Evaluation metrics
    tp <- length(intersect(pred_nonzero, nonzero_idx))
    fp <- length(setdiff(pred_nonzero, nonzero_idx))
    fn <- length(setdiff(nonzero_idx, pred_nonzero))
    tn <- p - tp - fp - fn
    cat("Run", i, "| TP:", tp, "FP:", fp, "FN:", fn, "TN:", tn, "\n")
    
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
    mcc_val <- mcc_custom(as.integer(true_beta != 0), as.integer(abs(beta) > 1e-4))
    rmse_val <- rmse(y, X_imp %*% beta)
    
    results[[i]] <- data.frame(Method = penalty, SupportSize = length(pred_nonzero),
                               F1 = f1, TPR = recall, FDR = fdr, MCC = mcc_val,
                               RMSE = rmse_val)
  }
  df <- do.call(rbind, results)
  total_runtime <- round(as.numeric(difftime(Sys.time(), overall_start, units = "secs")), 2)
  df$TotalRuntime <- total_runtime
  return(df)
}

#  Ensemble method (agreement ≥ 2) 
run_ensemble_real <- function(runs = 10) {
  cat(" Running Ensemble Method on Real Data\n")
  overall_start <- Sys.time()
  results <- vector("list", runs)
  
  for (i in 1:runs) {
    #cat("Run", i, "with ensemble\n")
    start_time <- Sys.time()
    
    # Step 1: Simulate beta and y
    true_beta <- rep(0, p)
    #nonzero_idx <- sample(1:p, 10)
    #true_beta[nonzero_idx] <- runif(10, -2, 2)
    nonzero_count <- min(10, p)
    nonzero_idx <- sample(1:p, nonzero_count)
    true_beta[nonzero_idx] <- runif(nonzero_count, -2, 2)
    
    y <- X %*% true_beta + rnorm(n)
    y <- scale(y)
    
    # Add outliers
    outlier_idx <- sample(1:n, floor(outlier_rate * n))
    y[outlier_idx] <- y[outlier_idx] + rnorm(length(outlier_idx), 0, contam_sd)
    
    # Step 2: Add missing values
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * missing_rate)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    
    # Step 3: Initial mean imputation
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    
    # Step 4: Ensemble via EM for SCAD, MCP, LOG
    penalties <- c("scad", "mcp", "log")
    beta_list <- list()
    for (penalty in penalties) {
      beta <- rep(0.01, p)
      for (em in 1:10) {
        X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
        colnames(X_df) <- paste0("V", 1:p)
        for (j in 1:p) {
          obs <- !is.na(X_miss[, j])
          if (sum(obs) < n) {
            fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))),
                          data = X_df[obs, ]), silent = TRUE)
            if (!inherits(fit, "try-error")) {
              X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
            }
          }
        }
        X_imp <- as.matrix(X_df)
        beta <- coordinate_descent(X_imp, y, beta, lambda, penalty)
      }
      beta_list[[penalty]] <- beta
    }
    
    # Step 5: Select features by agreement (≥2)
    agree_mask <- (abs(beta_list$scad) > 1e-4) + (abs(beta_list$mcp) > 1e-4) + (abs(beta_list$log) > 1e-4)
    selected_idx <- which(agree_mask >= 2)
    
    # Step 6: Final refit + evaluation
    beta_final <- rep(0, p)
    if (length(selected_idx) > 0) {
      fit <- try(lm(y ~ X_imp[, selected_idx] - 1), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        coefs <- coef(fit)
        if (any(is.na(coefs))) {
          coefs[is.na(coefs)] <- 0
        }
        beta_final[selected_idx] <- coefs
      }
    }
    
    pred <- X_imp %*% beta_final
    rmse_val <- if (any(is.na(pred))) NA else rmse(y, pred)
    
    # Metrics
    tp <- length(intersect(selected_idx, nonzero_idx))
    fp <- length(setdiff(selected_idx, nonzero_idx))
    fn <- length(setdiff(nonzero_idx, selected_idx))
    tn <- p - tp - fp - fn
    cat("Run", i, "| TP:", tp, "FP:", fp, "FN:", fn, "TN:", tn, "\n")
    
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
    
    denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc_val <- if (is.na(denom) || denom == 0) NA else (tp * tn - fp * fn) / denom
    
    #runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    results[[i]] <- data.frame(Method = "ensemble", SupportSize = length(selected_idx),
                               F1 = f1, TPR = recall, FDR = fdr, MCC = mcc_val,
                               RMSE = rmse_val)
  }
  
  # Combine and add total runtime
  df <- do.call(rbind, results)
  df$TotalRuntime <- round(as.numeric(difftime(Sys.time(), overall_start, units = "secs")), 2)
  return(df)
}


#  Run all methods 
results_all <- rbind(
  run_real_experiment("scad", runs),
  run_real_experiment("mcp", runs),
  run_real_experiment("log", runs),
  run_real_experiment("hybrid", runs),
  run_ensemble_real(runs)
)

#  Summarize: mean ± sd format 
summary_stats <- results_all %>%
  group_by(Method) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

# # add the toal time column:
# summary_formatted <- summary_stats %>%
#   mutate(
#     SupportSize = sprintf("%.1f ± %.1f", SupportSize_mean, SupportSize_sd),
#     F1          = sprintf("%.4f ± %.4f", F1_mean, F1_sd),
#     TPR         = sprintf("%.3f ± %.3f", TPR_mean, TPR_sd),
#     FDR         = sprintf("%.4f ± %.4f", FDR_mean, FDR_sd),
#     MCC         = sprintf("%.4f ± %.4f", MCC_mean, MCC_sd),
#     RMSE        = sprintf("%.2f ± %.2f", RMSE_mean, RMSE_sd)
#   ) %>%
#   select(Method, SupportSize, F1, TPR, FDR, MCC, RMSE) %>%
#   left_join(total_runtime, by = "Method") %>%
#   mutate(TotalTime = sprintf("%.2f sec", TotalTime))

### Not adding the total run time column:
summary_formatted <- summary_stats %>%
  mutate(
    SupportSize = sprintf("%.1f ± %.1f", SupportSize_mean, SupportSize_sd),
    F1          = sprintf("%.4f ± %.4f", F1_mean, F1_sd),
    TPR         = sprintf("%.3f ± %.3f", TPR_mean, TPR_sd),
    FDR         = sprintf("%.4f ± %.4f", FDR_mean, FDR_sd),
    MCC         = sprintf("%.4f ± %.4f", MCC_mean, MCC_sd),
    RMSE        = sprintf("%.2f ± %.2f", RMSE_mean, RMSE_sd)
  ) %>%
  select(Method, SupportSize, F1, TPR, FDR, MCC, RMSE)

#  save it
write.csv(summary_formatted, file = "summary_results.csv", row.names = FALSE)

######################### Real dataset 1 the top selected variables frequency under n runs:
# Load riboflavin data and prepare top 100 genes
library(hdi)
library(ggplot2)
data(riboflavin)
X_full <- scale(riboflavin$x)
top100 <- order(apply(X_full, 2, var), decreasing = TRUE)[1:100]
X <- X_full[, top100]
#gene_names <- colnames(riboflavin$x)[top100]
n <- nrow(X)
p <- ncol(X)
gene_names <- colnames(X)  
stopifnot(length(gene_names) == p) 

# Function to run real data and return selected variable indices
get_selected_vars <- function(penalty = "scad", runs = 100, lambda = 0.3) {
  selected_matrix <- matrix(0, nrow = runs, ncol = p)
  colnames(selected_matrix) <- gene_names
  
  
  for (i in 1:runs) {
    # Simulate sparse beta and y
    true_beta <- rep(0, p)
    nonzero_idx <- sample(1:p, min(10, p))
    true_beta[nonzero_idx] <- runif(length(nonzero_idx), -2, 2)
    y <- X %*% true_beta + rnorm(n)
    y <- as.numeric(scale(y))
    
    
    # Add missing and outliers
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * 0.2)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    y[sample(1:n, floor(0.1 * n))] <- y[sample(1:n, floor(0.1 * n))] + rnorm(floor(0.1 * n), 0, 25)
    
    # EM estimation
    beta <- rep(0.01, p)
    for (em in 1:10) {
      X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
      colnames(X_df) <- paste0("V", 1:p)
      for (j in 1:p) {
        obs <- !is.na(X_miss[, j])
        if (sum(obs) < n) {
          fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))), data = X_df[obs, ]), silent = TRUE)
          if (!inherits(fit, "try-error")) {
            X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
          }
        }
      }
      X_imp <- as.matrix(X_df)
      if (penalty == "hybrid") {
        beta <- coordinate_descent_hybrid(X_imp, y, beta, lambda)
      } else {
        beta <- coordinate_descent(X_imp, y, beta, lambda, penalty)
      }
    }
    
    selected_matrix[i, which(abs(beta) > 1e-4)] <- 1
  }
  
  return(colSums(selected_matrix))
}

# Function for ENSEMBLE method
get_selected_vars_ensemble <- function(runs = 100, lambda = 0.3) {
  freq <- numeric(p)
  names(freq) <- gene_names
  for (i in 1:runs) {
    true_beta <- rep(0, p)
    nonzero_idx <- sample(1:p, min(10, p))
    true_beta[nonzero_idx] <- runif(length(nonzero_idx), -2, 2)
    y <- X %*% true_beta + rnorm(n)
    y <- as.numeric(scale(y))
    
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * 0.2)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    y[sample(1:n, floor(0.1 * n))] <- y[sample(1:n, floor(0.1 * n))] + rnorm(floor(0.1 * n), 0, 25)
    
    penalties <- c("scad", "mcp", "log")
    beta_list <- list()
    for (penalty in penalties) {
      beta <- rep(0.01, p)
      for (em in 1:10) {
        X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
        colnames(X_df) <- paste0("V", 1:p)
        for (j in 1:p) {
          obs <- !is.na(X_miss[, j])
          if (sum(obs) < n) {
            fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))), data = X_df[obs, ]), silent = TRUE)
            if (!inherits(fit, "try-error")) {
              X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
            }
          }
        }
        X_imp <- as.matrix(X_df)
        beta <- coordinate_descent(X_imp, y, beta, lambda, penalty)
      }
      beta_list[[penalty]] <- beta
    }
    
    agree_mask <- (abs(beta_list$scad) > 1e-4) + (abs(beta_list$mcp) > 1e-4) + (abs(beta_list$log) > 1e-4)
    selected_idx <- which(agree_mask >= 2)
    freq[selected_idx] <- freq[selected_idx] + 1
  }
  return(freq)
}

# Get frequencies for all methods
freq_list <- list(
  SCAD     = get_selected_vars("scad", 100),
  MCP      = get_selected_vars("mcp", 100),
  LOG      = get_selected_vars("log", 100),
  HYBRID   = get_selected_vars("hybrid", 100),
  ENSEMBLE = get_selected_vars_ensemble(100)
)

# Plot for each method
library(ggplot2)

for (method in names(freq_list)) {
  df <- data.frame(Gene = names(freq_list[[method]]),
                   Frequency = freq_list[[method]])
  
  df <- df[order(-df$Frequency), ][1:10, ]  # Top 10 genes
  df$Gene <- factor(df$Gene, levels = df$Gene)  # Maintain order
  
  p <- ggplot(df, aes(x = Gene, y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = Frequency), vjust = -0.5, size = 4.5) +
    labs(x = NULL, y = "Frequency") +  # Removed title
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()  
    )
  
  print(p)
}

######################## real dataset 1: similarity for 10,20,30,40,50 runs:
get_selection_realdata <- function(runs = 10, lambda = 0.3, method = "scad") {
  selected_list <- list()
  for (i in 1:runs) {
    # simulate y
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, min(10, p))
    beta_true[nonzero_idx] <- runif(length(nonzero_idx), -2, 2)
    y <- X %*% beta_true + rnorm(n)
    y <- as.numeric(scale(y))
    
    # add missing and outliers
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * 0.2)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    y[sample(1:n, floor(0.1 * n))] <- y[sample(1:n, floor(0.1 * n))] + rnorm(floor(0.1 * n), 0, 25)
    
    # EM algorithm
    beta <- rep(0.01, p)
    for (em in 1:10) {
      X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
      colnames(X_df) <- paste0("V", 1:p)
      for (j in 1:p) {
        obs <- !is.na(X_miss[, j])
        if (sum(obs) < n) {
          fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))), data = X_df[obs, ]), silent = TRUE)
          if (!inherits(fit, "try-error")) {
            X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
          }
        }
      }
      X_imp <- as.matrix(X_df)
      if (method == "hybrid") {
        beta <- coordinate_descent_hybrid(X_imp, y, beta, lambda)
      } else {
        beta <- coordinate_descent(X_imp, y, beta, lambda, method)
      }
    }
    
    selected_list[[i]] <- which(abs(beta) > 1e-4)
  }
  return(selected_list)
}

get_ensemble_realdata <- function(runs = 10, lambda = 0.3) {
  selected_list <- list()
  for (i in 1:runs) {
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, min(10, p))
    beta_true[nonzero_idx] <- runif(length(nonzero_idx), -2, 2)
    y <- X %*% beta_true + rnorm(n)
    y <- as.numeric(scale(y))
    
    X_miss <- X
    miss_idx <- arrayInd(sample(1:(n * p), floor(n * p * 0.2)), .dim = c(n, p))
    X_miss[miss_idx] <- NA
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    y[sample(1:n, floor(0.1 * n))] <- y[sample(1:n, floor(0.1 * n))] + rnorm(floor(0.1 * n), 0, 25)
    
    penalties <- c("scad", "mcp", "log")
    beta_list <- list()
    for (penalty in penalties) {
      beta <- rep(0.01, p)
      for (em in 1:10) {
        X_df <- as.data.frame(matrix(unlist(X_imp), nrow = n, ncol = p))
        colnames(X_df) <- paste0("V", 1:p)
        for (j in 1:p) {
          obs <- !is.na(X_miss[, j])
          if (sum(obs) < n) {
            fit <- try(lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = "+"))), data = X_df[obs, ]), silent = TRUE)
            if (!inherits(fit, "try-error")) {
              X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
            }
          }
        }
        X_imp <- as.matrix(X_df)
        beta <- coordinate_descent(X_imp, y, beta, lambda, penalty)
      }
      beta_list[[penalty]] <- beta
    }
    
    agree_mask <- (abs(beta_list$scad) > 1e-4) + (abs(beta_list$mcp) > 1e-4) + (abs(beta_list$log) > 1e-4)
    selected_idx <- which(agree_mask >= 2)
    selected_list[[i]] <- selected_idx
  }
  return(selected_list)
}

# Run similarity analysis
run_similarity_analysis_realdata <- function(method = "scad", run_sizes = c(10, 20, 30, 40, 50), lambda = 0.3) {
  for (runs in run_sizes) {
    if (method == "ensemble") {
      selected_list <- get_ensemble_realdata(runs, lambda)
    } else {
      selected_list <- get_selection_realdata(runs, lambda, method)
    }
    
    similarities <- compute_jaccard_similarity(selected_list)
    cat("\nMethod:", toupper(method), " | Runs:", runs, "\n")
    cat("Mean:", mean(similarities), 
        "SD:", sd(similarities), 
        "Median:", median(similarities), 
        "Min:", min(similarities), 
        "Max:", max(similarities), "\n")
  }
}

# Example usage
run_similarity_analysis_realdata("scad")
run_similarity_analysis_realdata("mcp")
run_similarity_analysis_realdata("log")
run_similarity_analysis_realdata("hybrid")
run_similarity_analysis_realdata("ensemble")

##############################################
#Real data 2:breast Cancer Data
library(ggplot2)
library(dplyr)
library(tidyr)

#  Load the breast expression data 
X <- as.matrix(read.csv("breast_expression_original.csv"))
n <- nrow(X)
p <- ncol(X)

#  Set hyperparameters 
lambda <- 0.1  # adjust if needed
outlier_rate <- 0.05
contam_sd <- 5
missing_rate <- 0.05

#  Run all methods 
set.seed(123)
methods <- c("scad", "mcp", "log", "hybrid", "ensemble")
results_list <- list()

for (method in methods) {
  if (method == "ensemble") {
    results_list[[method]] <- run_ensemble_real(runs = 100)
  } else {
    results_list[[method]] <- run_real_experiment(penalty = method, runs = 100)
  }
}

combined_results <- do.call(rbind, lapply(names(results_list), function(m) {
  df <- results_list[[m]]
  df$Method <- m 
  df
}))


summary_table <- combined_results %>%
  group_by(Method) %>%
  summarise(
    SupportSize = sprintf("%.1f ± %.1f", mean(SupportSize), sd(SupportSize)),
    F1          = sprintf("%.4f ± %.4f", mean(F1), sd(F1)),
    TPR         = sprintf("%.3f ± %.3f", mean(TPR), sd(TPR)),
    FDR         = sprintf("%.4f ± %.4f", mean(FDR), sd(FDR)),
    MCC         = if (all(is.na(MCC))) "NA ± NA" else sprintf("%.4f ± %.4f", mean(MCC, na.rm = TRUE), sd(MCC, na.rm = TRUE)),
    RMSE        = sprintf("%.2f ± %.2f", mean(RMSE, na.rm = TRUE), sd(RMSE, na.rm = TRUE))
  )


write.csv(summary_table, "REAL_DATA_summary_results.csv", row.names = FALSE)

##########################plot the boxplot:
# Desired method order
method_order <- c("scad", "mcp", "log", "hybrid", "ensemble")

# Reshape data to long format
long_data <- combined_results %>%
  mutate(Method = factor(Method, levels = method_order)) %>%
  pivot_longer(cols = c(SupportSize, F1, TPR, FDR, MCC, RMSE),
               names_to = "Metric",
               values_to = "Value")

# Define metric names and filenames
metric_names <- unique(long_data$Metric)
file_labels <- letters[1:length(metric_names)]  # a, b, c, d, e, f

# Loop and save each metric as pdf
for (i in seq_along(metric_names)) {
  metric_name <- metric_names[i]
  file_name <- paste0("Fig1", file_labels[i], ".pdf")
  metric_data <- long_data %>% filter(Metric == metric_name)
  
  p <- ggplot(metric_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      legend.position = "none",
      legend.text = element_text(size = 20),                         # Prepare legend font size
      legend.title = element_text(size = 20)
    ) +
    labs(x = "Method", y = metric_name)
  
  ggsave(filename = file_name, plot = p, width = 7, height = 5)
}

############ real dataset 2 similarity:
# Breast Cancer Dataset Similarity Analysis
library(dplyr)

# Load Breast Cancer gene expression data
X <- as.matrix(read.csv("breast_expression_original.csv"))
n <- nrow(X)
p <- ncol(X)

# parameters
lambda <- 0.3
run_sizes <- c(10, 20, 30, 40, 50)
set.seed(123)

# we still use the same helper funcs of real data 1 here.

# Run All Methods
methods <- c("scad", "mcp", "log", "hybrid", "ensemble")
for (m in methods) {
  run_similarity_analysis_realdata(method = m, run_sizes = run_sizes, lambda = lambda)
}
