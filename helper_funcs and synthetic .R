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
    
    # Correct split (handles Total_Runtime)
    res_split <- res %>%
      mutate(
        Mean = ifelse(Metric == "Total_Runtime",
                      as.numeric(Mean_SD),
                      as.numeric(sub(" ±.*", "", Mean_SD))),
        SD   = ifelse(Metric == "Total_Runtime",
                      0,
                      as.numeric(sub(".*± ", "", Mean_SD))),
        Lambda = lam,
        Penalty = penalty
      ) %>%
      select(Metric, Mean, SD, Lambda, Penalty)
    
    lambda_results[[paste0(penalty, "_", lam)]] <- res_split
  }
}

df_lambda <- do.call(rbind, lambda_results)
saveRDS(df_lambda, "df_lambda_results.rds")

df_lambda_clean <- df_lambda  # <-- no further mutate needed!
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
# penalty_linetypes <- c(
#   "scad"     = "dashed",
#   "mcp"      = "dashed",
#   "log"      = "dashed",
#   "hybrid"   = "dashed",
#   "ensemble" = "dashed"
# )                         
penalty_linetypes <- c(
  "scad"     = "solid",
  "mcp"      = "longdash",
  "log"      = "twodash",
  "hybrid"   = "dotdash",
  "ensemble" = "dotted"
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

#####################################
##############Below Revision 11/17/2025:
# Additional Noise Experiments (Reviewer 2 #6)
#Some clarification of the large RMSE values reported for 
# SCAD, MCP, LOG, and HYBRID in Table 2. To address this, we run an additional 
# diagnostic experiment under the same simulation setting as Table 2. 
# The goal is to illustrate how coefficient variability affects RMSE across 
# penalties and to confirm that ENSEMBLE maintains much smaller and stable RMSE.
#
# Experiment Setup:
#   • n = 100 observations, p = 100 predictors
#   • s = 10 true nonzero coefficients
#   • X generated from AR(1) covariance (rho = 0.5)
#   • 20% missing entries in X (MCAR), imputed in EM steps
#   • 10% responses contaminated by N(0, 25^2)
#   • Penalties evaluated: SCAD, MCP, LOG, HYBRID, ENSEMBLE
#   • RMSE computed on the s = 10 true active coefficients
#   • Repeated for 3 runs (sufficient for diagnostic summary)
#
# Output:
# The summary table created by this script is used as Table 3 in the 
# manuscript revision. It shows that single-penalty methods yield RMSE 
# around 300 due to noisy false positives, while ENSEMBLE achieves RMSE 
# around 2–3 due to consensus-based variable selection.

generate_noise <- function(type, n) {
  if (type == "t6") {
    return(rt(n, df = 6))
  }
  if (type == "gauss5") {
    base <- rnorm(n)
    contam <- rbinom(n, size = 1, prob = 0.05)
    return(base + contam * rnorm(n, 0, 10))
  }
  if (type == "heavy") {
    return(rnorm(n, 0, 25))
  }
}

run_experiment_noise <- function(noise_type = "t6", runs = 3, lambda = 0.3) {
  
  cat("\n")
  cat("Running noise scenario:", noise_type, "\n")
  
  results_scad <- list()
  results_mcp  <- list()
  results_log  <- list()
  results_hybrid <- list()
  results_ensemble <- list()
  
  for (run in 1:runs) {
    
    X_full <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    beta_true <- rep(0, p)
    nonzero_idx <- sample(1:p, s)
    beta_true[nonzero_idx] <- runif(s, 1, 2) * sample(c(-1, 1), s, TRUE)
    
    eps <- generate_noise(noise_type, n)
    y <- X_full %*% beta_true + eps
    
    X_miss <- X_full
    num_missing <- floor(missing_rate * n * p)
    miss_idx <- arrayInd(sample(1:(n*p), num_missing), .dim = c(n, p))
    for (k in 1:nrow(miss_idx)) {
      X_miss[miss_idx[k, 1], miss_idx[k, 2]] <- NA
    }
    
    results_scad[[run]] <- as.numeric(
      run_em_algorithm(
        X_miss, y, beta_true, nonzero_idx,
        lambda = lambda, penalty = "scad"
      )[1:5]
    )
    
    results_mcp[[run]] <- as.numeric(
      run_em_algorithm(
        X_miss, y, beta_true, nonzero_idx,
        lambda = lambda, penalty = "mcp"
      )[1:5]
    )
    
    results_log[[run]] <- as.numeric(
      run_em_algorithm(
        X_miss, y, beta_true, nonzero_idx,
        lambda = lambda, penalty = "log"
      )[1:5]
    )
    
    results_hybrid[[run]] <- as.numeric(
      run_em_hybrid(
        X_miss, y, beta_true, nonzero_idx,
        lambda = lambda
      )[1:5]
    )
    
    ens <- run_ensemble_experiment(lambda = lambda, runs = 1)
    
    metrics_needed <- c("F1", "TPR", "FDR", "MCC", "RMSE")
    ens_vals <- sapply(metrics_needed, function(m) {
      val <- ens$Mean_SD[ens$Metric == m]
      as.numeric(sub(" ±.*", "", val))
    })
    
    results_ensemble[[run]] <- as.numeric(ens_vals)
  }
  
  df_scad <- do.call(rbind, results_scad)
  df_mcp <- do.call(rbind, results_mcp)
  df_log <- do.call(rbind, results_log)
  df_hybrid <- do.call(rbind, results_hybrid)
  df_ensemble <- do.call(rbind, results_ensemble)
  
  summary_table <- rbind(
    cbind(Method="SCAD",     
          setNames(t(colMeans(df_scad)), 
                   c("F1","TPR","FDR","MCC","RMSE"))),
    
    cbind(Method="MCP",      
          setNames(t(colMeans(df_mcp)),  
                   c("F1","TPR","FDR","MCC","RMSE"))),
    
    cbind(Method="LOG",      
          setNames(t(colMeans(df_log)),  
                   c("F1","TPR","FDR","MCC","RMSE"))),
    
    cbind(Method="HYBRID",   
          setNames(t(colMeans(df_hybrid)), 
                   c("F1","TPR","FDR","MCC","RMSE"))),
    
    cbind(Method="ENSEMBLE",
          setNames(t(colMeans(df_ensemble)), 
                   c("F1","TPR","FDR","MCC","RMSE")))
  )
  
  print(summary_table)
  return(summary_table)
}
#run it
summary_table_t6 <- run_experiment_noise("t6", runs = 3)
summary_table_gauss5 <- run_experiment_noise("gauss5", runs = 3)
summary_table_heavy <- run_experiment_noise("heavy", runs = 3)

#write it                        
write.csv(summary_table_t6, "NoiseResults_t6.csv", row.names = FALSE)
write.csv(summary_table_gauss5, "NoiseResults_gauss5.csv", row.names = FALSE)
write.csv(summary_table_heavy, "NoiseResults_heavy.csv", row.names = FALSE)
                         
