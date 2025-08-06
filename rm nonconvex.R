#variable selection paper:
#################################
#Synthetic dataset:
##################################
#Using 3 penalties(SCAD,MCP,LOG) get 6 metrics:
# RM-Nonconvex with SCAD, MCP, and Log Penalties
library(MASS)
library(Matrix)
library(openxlsx)
library(Metrics)

set.seed(123)
n <- 100; p <- 100; s <- 10; rho <- 0.5
missing_rate <- 0.2; outlier_rate <- 0.1
sigma <- 1; contam_sd <- 25

Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))

# Penalty Derivatives
scad_deriv <- function(beta, lambda, a = 3.7) {
  absb <- abs(beta)
  out <- ifelse(absb <= lambda, lambda,
                ifelse(absb <= a * lambda, (a * lambda - absb) / (a - 1), 0))
  out * sign(beta)
}

mcp_deriv <- function(beta, lambda, gamma = 3) {
  abs_beta <- abs(beta)
  out <- ifelse(abs_beta <= gamma * lambda,
                lambda * (1 - abs_beta / (gamma * lambda)),
                0)
  return(out * sign(beta))
}

log_deriv <- function(beta, lambda, epsilon = 0.1) {
  return(lambda / (epsilon + abs(beta)) * sign(beta))
}

mcc_custom <- function(y_true, y_pred) {
  TP <- sum(y_true == 1 & y_pred == 1)
  TN <- sum(y_true == 0 & y_pred == 0)
  FP <- sum(y_true == 0 & y_pred == 1)
  FN <- sum(y_true == 1 & y_pred == 0)
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (is.na(denom) || denom == 0) return(0)
  (TP * TN - FP * FN) / denom
}

coordinate_descent <- function(X, y, beta_init, lambda, penalty = "scad", delta = 1.345, max_iter = 50, step_size = 0.001, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- beta_init
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    
    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      z_j <- X[, j]
      r_j_huber <- ifelse(abs(r_j) <= delta, r_j, delta * sign(r_j))
      grad <- -sum(z_j * r_j_huber) / n
      
      # Penalty derivative
      if (penalty == "scad") {
        w_j <- scad_deriv(beta[j], lambda)
      } else if (penalty == "mcp") {
        w_j <- mcp_deriv(beta[j], lambda)
      } else if (penalty == "log") {
        w_j <- log_deriv(beta[j], lambda)
      } else {
        stop("Unknown penalty")
      }
      
      grad <- max(min(grad, 10), -10)
      beta[j] <- beta[j] - step_size * (grad + w_j)
      beta[j] <- sign(beta[j]) * min(abs(beta[j]), 100)
    }
    
    if (sqrt(sum((beta - beta_old)^2)) < tol) break
  }
  
  return(beta)
  
}


run_em_algorithm <- function(X_miss, y, beta_true, nonzero_idx,
                             lambda = 0.3, max_em = 10, penalty = "SCAD") {
  n <- nrow(X_miss); p <- ncol(X_miss)
  X_imp <- X_miss
  for (j in 1:p) {
    X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
  }
  beta <- rep(0.01, p)
  start_time <- Sys.time()
  for (em in 1:max_em) {
    X_df <- as.data.frame(X_imp); colnames(X_df) <- paste0("V", 1:p)
    for (j in 1:p) {
      obs <- !is.na(X_miss[, j])
      if (sum(obs) < n) {
        fit <- lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = " + "))),
                  data = X_df[obs, ])
        X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
      }
    }
    X_imp <- as.matrix(X_df)
    beta <- coordinate_descent(X_imp, y, beta, lambda, penalty = penalty)
  }
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  pred_nonzero <- which(abs(beta) > 1e-4)
  tp <- length(intersect(pred_nonzero, nonzero_idx))
  fp <- length(setdiff(pred_nonzero, nonzero_idx))
  fn <- length(setdiff(nonzero_idx, pred_nonzero))
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  mcc_val <- mcc_custom(as.integer(beta_true != 0), as.integer(abs(beta) > 1e-4))
  rmse_val <- rmse(y, X_imp %*% beta)
  data.frame(F1 = f1, TPR = recall, FDR = fdr, MCC = mcc_val, RMSE = rmse_val, Runtime = runtime)
}

#  Wrapper function to switch penalty
run_experiment_with_penalty <- function(penalty_type, lambda = 0.3, runs = 50) {
  cat("\nRunning with", penalty_type, "penalty\n")
  overall_start <- Sys.time()
  results <- vector("list", runs)
  
  for (run in 1:runs) {
    #cat("Run", run, "...\n")  not plan to catch the running clock here
    
    # regenerate data each time
    X_full <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2) * sample(c(-1, 1), s, replace = TRUE)
    y_clean <- X_full %*% beta_true + rnorm(n, 0, sigma)
    y <- y_clean
    outlier_idx <- sample(1:n, floor(outlier_rate * n))
    y[outlier_idx] <- y[outlier_idx] + rnorm(length(outlier_idx), 0, contam_sd)
    
    X_miss <- X_full
    num_missing <- floor(missing_rate * n * p)
    miss_idx <- arrayInd(sample(1:(n * p), num_missing), .dim = c(n, p))
    for (k in 1:nrow(miss_idx)) {
      X_miss[miss_idx[k, 1], miss_idx[k, 2]] <- NA
    }
    
    # run method with selected penalty
    results[[run]] <- run_em_algorithm(X_miss, y, beta_true, nonzero_idx,
                                       lambda = lambda, penalty = penalty_type)
  }
  
  df <- do.call(rbind, results)
  total_runtime <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
  
  # Compute mean and SD for the first 5 metrics
  avg <- colMeans(df[, 1:5], na.rm = TRUE)
  sds <- apply(df[, 1:5], 2, sd, na.rm = TRUE)
  
  # Format mean ± sd
  summary_formatted <- paste0(
    round(avg, 4), " ± ", round(sds, 4)
  )
  
  # Build summary table
  summary_stats <- data.frame(
    Metric = names(avg),
    `Mean ± SD` = summary_formatted
  )
  
  # Add Total Runtime as a separate row
  summary_stats <- rbind(
    summary_stats,
    data.frame(Metric = "Total_Runtime", `Mean ± SD` = round(total_runtime, 4))
  )
  
  cat("\nSummary (Mean ± SD) for", penalty_type, "\n")
  print(summary_stats)
  return(summary_stats)
  
}


#  Run experiments for each penalty :
results_scad <- run_experiment_with_penalty("scad",runs = 50)
results_mcp  <- run_experiment_with_penalty("mcp",runs = 50)
results_log  <- run_experiment_with_penalty("log",runs = 50)

################The hybrid method implement:
# Hybrid Penalty Derivative 
hybrid_deriv <- function(beta, lambda, alpha = c(1/3, 1/3, 1/3), a = 3.7, gamma = 3, epsilon = 0.1) {
  # SCAD
  absb <- abs(beta)
  scad <- ifelse(absb <= lambda, lambda,
                 ifelse(absb <= a * lambda, (a * lambda - absb) / (a - 1), 0))
  scad <- scad * sign(beta)
  
  # MCP
  mcp <- ifelse(absb <= gamma * lambda,
                lambda * (1 - absb / (gamma * lambda)),
                0)
  mcp <- mcp * sign(beta)
  
  # Log
  logpen <- lambda / (epsilon + absb) * sign(beta)
  
  # Combine
  alpha <- alpha / sum(alpha)  # normalize
  hybrid <- alpha[1] * scad + alpha[2] * mcp + alpha[3] * logpen
  return(hybrid)
}

#  Coordinate Descent with Hybrid Penalty 
coordinate_descent_hybrid <- function(X, y, beta_init, lambda, alpha = c(1/3, 1/3, 1/3),
                                      delta = 1.345, max_iter = 50, step_size = 0.001, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- beta_init
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      z_j <- X[, j]
      r_j_huber <- ifelse(abs(r_j) <= delta, r_j, delta * sign(r_j))
      grad <- -sum(z_j * r_j_huber) / n
      w_j <- hybrid_deriv(beta[j], lambda, alpha)
      grad <- max(min(grad, 10), -10)
      beta[j] <- beta[j] - step_size * (grad + w_j)
      beta[j] <- sign(beta[j]) * min(abs(beta[j]), 100)
    }
    if (sqrt(sum((beta - beta_old)^2)) < tol) break
  }
  return(beta)
  
}

#  EM Algorithm with Hybrid Penalty 
run_em_hybrid <- function(X_miss, y, beta_true, nonzero_idx,
                          lambda = 0.3, alpha = c(1/3, 1/3, 1/3), max_em = 10) {
  n <- nrow(X_miss); p <- ncol(X_miss)
  X_imp <- X_miss
  for (j in 1:p) {
    X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
  }
  beta <- rep(0.01, p)
  start_time <- Sys.time()
  for (em in 1:max_em) {
    X_df <- as.data.frame(X_imp); colnames(X_df) <- paste0("V", 1:p)
    for (j in 1:p) {
      obs <- !is.na(X_miss[, j])
      if (sum(obs) < n) {
        fit <- lm(as.formula(paste0("V", j, " ~ ", paste0("V", setdiff(1:p, j), collapse = " + "))),
                  data = X_df[obs, ])
        X_df[!obs, j] <- predict(fit, newdata = X_df[!obs, ])
      }
    }
    X_imp <- as.matrix(X_df)
    beta <- coordinate_descent_hybrid(X_imp, y, beta, lambda, alpha)
  }
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (!is.numeric(beta)) beta <- as.numeric(unlist(beta))
  beta[is.na(beta)] <- 0
  pred_nonzero <- which(abs(beta) > 1e-4)
  tp <- length(intersect(pred_nonzero, nonzero_idx))
  fp <- length(setdiff(pred_nonzero, nonzero_idx))
  fn <- length(setdiff(nonzero_idx, pred_nonzero))
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  mcc_val <- mcc_custom(as.integer(beta_true != 0), as.integer(abs(beta) > 1e-4))
  rmse_val <- rmse(y, X_imp %*% beta)
  data.frame(F1 = f1, TPR = recall, FDR = fdr, MCC = mcc_val, RMSE = rmse_val, Runtime = runtime)
}

#  Run Hybrid Experiment:
run_hybrid_experiment <- function(lambda = 0.3, alpha = c(1/3, 1/3, 1/3), runs = 50) {
  cat("\n Running Hybrid Method \n")
  overall_start <- Sys.time()  # Start total timing
  results <- vector("list", runs)
  
  for (run in 1:runs) {
    #cat("Run", run, "...\n")
    X_full <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2) * sample(c(-1, 1), s, replace = TRUE)
    y_clean <- X_full %*% beta_true + rnorm(n, 0, sigma)
    y <- y_clean
    outlier_idx <- sample(1:n, floor(outlier_rate * n))
    y[outlier_idx] <- y[outlier_idx] + rnorm(length(outlier_idx), 0, contam_sd)
    X_miss <- X_full
    num_missing <- floor(missing_rate * n * p)
    miss_idx <- arrayInd(sample(1:(n * p), num_missing), .dim = c(n, p))
    for (k in 1:nrow(miss_idx)) {
      X_miss[miss_idx[k, 1], miss_idx[k, 2]] <- NA
    }
    results[[run]] <- run_em_hybrid(X_miss, y, beta_true, nonzero_idx,
                                    lambda = lambda, alpha = alpha)
  }
  
  overall_end <- Sys.time()
  total_runtime <- as.numeric(difftime(overall_end, overall_start, units = "secs"))
  
  df <- do.call(rbind, results)
  means <- colMeans(df[, 1:5], na.rm = TRUE)
  sds   <- apply(df[, 1:5], 2, sd, na.rm = TRUE)
  
  summary_df <- data.frame(
    Metric = names(means),
    Mean_SD = sprintf("%.4f ± %.4f", means, sds)
  )
  
  # Add total runtime as a separate row
  summary_df <- rbind(summary_df, data.frame(Metric = "Total_Runtime", 
                                             Mean_SD = sprintf("%.3f", total_runtime)))
  
  cat("\nSummary (Mean ± SD) for Hybrid Method\n")
  print(summary_df, row.names = FALSE)
  
  return(summary_df)
}
results_hybrid <- run_hybrid_experiment(runs = 50)

#Run Ensember Experiment:
run_ensemble_experiment <- function(lambda = 0.3, runs = 50) {
  cat(" Running Ensemble Method (Agreement of ≥2) \n")
  overall_start <- Sys.time()  # Total timing start
  results <- vector("list", runs)
  
  for (run in 1:runs) {
    #cat("Run", run, "...\n")
    start_time <- Sys.time()
    
    # Data generation 
    X_full <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2) * sample(c(-1, 1), s, replace = TRUE)
    y_clean <- X_full %*% beta_true + rnorm(n, 0, sigma)
    y <- y_clean
    outlier_idx <- sample(1:n, floor(outlier_rate * n))
    y[outlier_idx] <- y[outlier_idx] + rnorm(length(outlier_idx), 0, contam_sd)
    
    X_miss <- X_full
    num_missing <- floor(missing_rate * n * p)
    miss_idx <- arrayInd(sample(1:(n * p), num_missing), .dim = c(n, p))
    for (k in 1:nrow(miss_idx)) {
      X_miss[miss_idx[k, 1], miss_idx[k, 2]] <- NA
    }
    
    # Impute 
    X_imp <- X_miss
    for (j in 1:p) {
      X_imp[is.na(X_imp[, j]), j] <- mean(X_imp[, j], na.rm = TRUE)
    }
    
    #  Coordinate descent 
    beta_scad <- try(coordinate_descent(X_imp, y, rep(0.01, p), lambda, penalty = "scad"), silent = TRUE)
    beta_mcp  <- try(coordinate_descent(X_imp, y, rep(0.01, p), lambda, penalty = "mcp"),  silent = TRUE)
    beta_log  <- try(coordinate_descent(X_imp, y, rep(0.01, p), lambda, penalty = "log"),  silent = TRUE)
    
    if (inherits(beta_scad, "try-error") || inherits(beta_mcp, "try-error") || inherits(beta_log, "try-error")) {
      warning("Skipping run due to error.")
      next
    }
    
    beta_scad_vec <- abs(as.numeric(beta_scad))
    beta_mcp_vec  <- abs(as.numeric(beta_mcp))
    beta_log_vec  <- abs(as.numeric(beta_log))
    
    select_mask <- (beta_scad_vec > 1e-4) + (beta_mcp_vec > 1e-4) + (beta_log_vec > 1e-4) >= 2
    selected_idx <- which(select_mask)
    
    if (length(selected_idx) == 0) {
      warning("No features selected in run ", run)
      next
    }
    
    beta_final <- rep(0, p)
    fit <- lm(y ~ X_imp[, selected_idx] - 1)
    beta_final[selected_idx] <- coef(fit)
    
    tp <- length(intersect(selected_idx, nonzero_idx))
    fp <- length(setdiff(selected_idx, nonzero_idx))
    fn <- length(setdiff(nonzero_idx, selected_idx))
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
    mcc_val <- mcc_custom(as.integer(beta_true != 0), as.integer(abs(beta_final) > 1e-4))
    rmse_val <- tryCatch({
      rmse(y, X_imp %*% beta_final)
    }, error = function(e) NA)
    
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results[[run]] <- data.frame(F1 = f1, TPR = recall, FDR = fdr, 
                                 MCC = mcc_val, RMSE = rmse_val, Runtime = runtime)
  }
  
  overall_end <- Sys.time()
  total_runtime <- as.numeric(difftime(overall_end, overall_start, units = "secs"))
  
  df <- do.call(rbind, results)
  means <- colMeans(df, na.rm = TRUE)
  sds <- apply(df, 2, sd, na.rm = TRUE)
  
  summary_df <- data.frame(
    Metric = names(means),
    Mean_SD = sprintf("%.4f ± %.4f", means, sds)
  )
  
  # Add total runtime as final row
  summary_df <- rbind(summary_df, data.frame(Metric = "Total_Runtime", 
                                             Mean_SD = sprintf("%.3f", total_runtime)))
  
  cat("\nSummary (Mean ± SD) for Ensemble Method\n")
  print(summary_df, row.names = FALSE)
  
  return(summary_df)
}

results_ensemble <- run_ensemble_experiment(runs = 50)

#################### Different Lambda Curves for Each Penalty(synthetic data) 
# I used runs=10 here
lambda_seq <- c(0.05, 0.3, 0.6, 0.9)
penalties <- c("scad", "mcp", "log", "hybrid", "ensemble")
lambda_results <- list()

for (penalty in penalties) {
  for (lam in lambda_seq) {
    cat("\nRunning", penalty, "with lambda =", lam, "\n")
    
    if (penalty %in% c("scad", "mcp", "log")) {
      res <- run_experiment_with_penalty(penalty, lambda = lam, runs = 10)
    } else if (penalty == "hybrid") {
      res <- run_hybrid_experiment(lambda = lam, runs = 10)
    } else if (penalty == "ensemble") {
      res <- run_ensemble_experiment(lambda = lam, runs = 10)
    }
    
    if ("Mean...SD" %in% colnames(res)) {
      names(res)[names(res) == "Mean...SD"] <- "Mean_SD"
    }
    
    # Split Mean ± SD into separate columns
    res_split <- res %>%
      mutate(
        Mean = as.numeric(sub(" ±.*", "", Mean_SD)),
        SD   = as.numeric(sub(".*± ", "", Mean_SD)),
        Lambda = lam,
        Penalty = penalty
      ) %>%
      select(Metric, Mean, SD, Lambda, Penalty)
    
    lambda_results[[paste0(penalty, "_", lam)]] <- res_split
  }
}

# Combine and save
df_lambda <- do.call(rbind, lambda_results)
saveRDS(df_lambda, "df_lambda_results.rds")

df_lambda_clean <- df_lambda %>%
  filter(Metric != "Total_Runtime")

unique_metrics <- unique(df_lambda_clean$Metric)
penalty_levels <- c("scad", "mcp", "log", "hybrid", "ensemble")

# Labels shown in the legend
penalty_labels <- c(
  "scad"     = "SCAD",
  "mcp"      = "MCP",
  "log"      = "LOG",
  "hybrid"   = "HYBRID",
  "ensemble" = "ENSEMBLE"
)

# Colors for each penalty
penalty_colors <- c(
  "scad"     = "red",
  "mcp"      = "orange",
  "log"      = "blue",
  "hybrid"   = "black",
  "ensemble" = "brown"
)

# All dashed line styles
penalty_linetypes <- c(
  "scad"     = "dashed",
  "mcp"      = "dashed",
  "log"      = "dashed",
  "hybrid"   = "dashed",
  "ensemble" = "dashed"
)


for (metric in unique_metrics) {
  plot_data <- df_lambda_clean %>% filter(Metric == metric)
  plot_data$Penalty <- factor(plot_data$Penalty, levels = penalty_levels)
  
  p <- ggplot(plot_data, aes(x = Lambda, y = Mean, color = Penalty, linetype = Penalty)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = penalty_colors, labels = penalty_labels) +
    scale_linetype_manual(values = penalty_linetypes, labels = penalty_labels) +
    labs(
      x = expression(lambda),
      y = metric
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = "top",
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20)
    )
  
  ggsave(filename = paste0("LambdaPlot_", metric, ".pdf"), plot = p, width = 7, height = 5)
}

########################## Synthetic data similarities of adjacent runs:
# Jaccard Similarity Function 
compute_jaccard_similarity <- function(selection_list) {
  n <- length(selection_list)
  similarities <- numeric(n - 1)
  for (i in 1:(n - 1)) {
    set1 <- selection_list[[i]]
    set2 <- selection_list[[i + 1]]
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    similarities[i] <- if (union == 0) 1 else intersection / union
  }
  return(similarities)
}

#  Penalty-Specific Selection 
get_selection <- function(runs = 10, lambda = 0.3, method = "scad") {
  selected_list <- list()
  for (run in 1:runs) {
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2)
    y <- X %*% beta_true + rnorm(n)
    X_imp <- X
    beta_init <- rep(0.01, p)
    
    beta_hat <- coordinate_descent(X_imp, y, beta_init, lambda, penalty = method)
    selected_list[[run]] <- which(abs(beta_hat) > 1e-4)
  }
  return(selected_list)
}

#  Hybrid Selection 
get_hybrid_selection <- function(runs = 10, lambda = 0.3) {
  selected_list <- list()
  for (run in 1:runs) {
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2)
    y <- X %*% beta_true + rnorm(n)
    X_imp <- X
    beta_init <- rep(0.01, p)
    beta_hat <- coordinate_descent_hybrid(X_imp, y, beta_init, lambda)
    selected_list[[run]] <- which(abs(beta_hat) > 1e-4)
  }
  return(selected_list)
}

# Ensemble Selection 
get_ensemble_selection <- function(runs = 10, lambda = 0.3) {
  selected_list <- list()
  for (run in 1:runs) {
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2)
    y <- X %*% beta_true + rnorm(n)
    X_imp <- X
    beta_init <- rep(0.01, p)
    
    beta_scad <- coordinate_descent(X_imp, y, beta_init, lambda, penalty = "scad")
    beta_mcp  <- coordinate_descent(X_imp, y, beta_init, lambda, penalty = "mcp")
    beta_log  <- coordinate_descent(X_imp, y, beta_init, lambda, penalty = "log")
    
    select_mask <- (abs(beta_scad) > 1e-4) + (abs(beta_mcp) > 1e-4) + (abs(beta_log) > 1e-4) >= 2
    selected_idx <- which(select_mask)
    selected_list[[run]] <- selected_idx
  }
  return(selected_list)
}

# Run Analysis 
run_similarity_analysis <- function(method = "scad", run_sizes = c(10, 20, 30, 40, 50), lambda = 0.3) {
  cat("\nRunning for penalty:", method, "\n")
  for (runs in run_sizes) {
    if (method == "hybrid") {
      selected_list <- get_hybrid_selection(runs, lambda)
    } else if (method == "ensemble") {
      selected_list <- get_ensemble_selection(runs, lambda)
    } else {
      selected_list <- get_selection(runs, lambda, method)
    }
    
    similarities <- compute_jaccard_similarity(selected_list)
    summary_stats <- data.frame(
      Runs = runs,
      Mean = mean(similarities),
      SD = sd(similarities),
      Median = median(similarities),
      Min = min(similarities),
      Max = max(similarities)
    )
    print(summary_stats)
  }
}

run_similarity_analysis("scad")
run_similarity_analysis("mcp")
run_similarity_analysis("log")
run_similarity_analysis("hybrid")
run_similarity_analysis("ensemble")

#######################Synthetic data ablation study:
run_ablation_hybrid_experiment <- function(lambda = 0.3, runs = 50) {
  ablation_settings <- list(
    "SCAD + MCP + LOG" = c(1/3, 1/3, 1/3),
    "SCAD + MCP"       = c(0.5, 0.5, 0),
    "SCAD + LOG"       = c(0.5, 0, 0.5),
    "MCP + LOG"        = c(0, 0.5, 0.5)
  )
  
  ablation_results <- list()
  
  for (label in names(ablation_settings)) {
    alpha <- ablation_settings[[label]]
    cat("\nRunning Ablation:", label, "\n")
    
    result <- run_hybrid_experiment(lambda = lambda, alpha = alpha, runs = runs)
    result$Setting <- label
    ablation_results[[label]] <- result
  }
  
  # Combine into a single summary table
  df_ablation <- do.call(rbind, ablation_results)
  df_ablation <- df_ablation[, c("Setting", "Metric", "Mean_SD")]
  df_ablation <- reshape(df_ablation, timevar = "Metric", idvar = "Setting", direction = "wide")
  colnames(df_ablation) <- gsub("Mean_SD\\.", "", colnames(df_ablation))
  return(df_ablation)
}

# Run and view results
ablation_summary <- run_ablation_hybrid_experiment(runs = 50)
print(ablation_summary)


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

######################### Real dataset 1 the top selected variables frequency under n runs(I didn't put it in the paper):
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
