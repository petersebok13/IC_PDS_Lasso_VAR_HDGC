TM_check_metrics <- function(est_beta, true_beta, tol = 0.0001) {
  
  est_beta[is.na(est_beta)] <- 0     # Treat NAs as zero
  true_beta[is.na(true_beta)] <- 0   # Very safe
  est_beta[abs(est_beta) < tol] <- 0
  if (length(est_beta) != length(true_beta)) {
    stop("Mismatch between estimated and true beta dimensions.")
  }
  
  
  S     <- which(true_beta != 0)     # True driver variables
  S_c   <- which(true_beta == 0)     # True non-drivers
  S_est <- which(est_beta != 0)      # Selected variables
  
  TD <- intersect(S, S_est)          # True discoveries
  FD <- intersect(S_c, S_est)        # False discoveries
  
  
  accuracy <- TRUE
  for (i in 1:length(true_beta)) {
    if (true_beta[i] != 0) {
      if (est_beta[i] == 0) {
        accuracy = FALSE
      }
    } else {
      if (est_beta[i] != 0) {
        accuracy = FALSE
      }
    }
  }
  FDP      <- length(FD) / max(length(S_est), 1)  # False Discovery Proportion
  TPP      <- length(TD) / max(length(S), 1)      # True Positive Proportion
  
  U_fit = ifelse(TPP < 1, TRUE, FALSE)
  O_fit = ifelse(TPP == 1 && FDP > 0, TRUE, FALSE)
  
  return(list(accuracy = accuracy, U_fit = U_fit, O_fit = O_fit))
}


################################################################################
TM_score_test_VAR <- function(Beta, n, ntest, lag = 1, varcov = "Identity") {
  
  # Initialize metrics
  methods <- c("Theo", "CV", "AIC", "BIC", "EBIC", "ERIC")
  acc <- setNames(rep(0, length(methods)), methods)
  U_fit <- setNames(rep(0, length(methods)), methods)
  O_fit <- setNames(rep(0, length(methods)), methods)
  
  for (i in 1:ntest) {
    data <- gen_sim_data_ts(Beta_matrix = Beta, n = n, varcov = varcov, lag = 1)
    
    # --- Theoretical lambda
    sigma_squared_OLS <- mean((lm(data$y ~ data$x)$residual)^2)
    lambda_theo <- sqrt(sigma_squared_OLS) / sqrt(n) * qnorm(1 - 0.05 / (2 * ncol(data$x)*log(nrow(data$x)))) #check
    coef_theo <- glmnet(data$x, data$y, alpha = 1, lambda = lambda_theo)
    coef_theo <- as.numeric(coef(coef_theo))[-1]
    coef_theo[is.na(coef_theo)] <- 0
    met_theo <- TM_check_metrics(coef_theo, Beta[1,])
    
    # --- Cross-Validation
    coef_cv <- rolling_cv_glmnet(data$x, data$y, alpha = 1, train_ratio = 0.8)
    coef_cv <- as.numeric(coef_cv$final_model$beta)
    coef_cv[is.na(coef_cv)] <- 0
    met_cv <- TM_check_metrics(t(coef_cv), Beta[1,])
    
    # --- AIC
    coef_aic <- ic_glmnetboundT(data$x, data$y, crit = "aic")$coefficients[-1]
    met_aic <- TM_check_metrics(coef_aic, Beta[1,])
    
    # --- BIC
    coef_bic <- ic_glmnetboundT(data$x, data$y, crit = "bic")$coefficients[-1]
    met_bic <- TM_check_metrics(coef_bic, Beta[1,])
    
    # --- EBIC
    coef_ebic <- ic_glmnetboundT(data$x, data$y, crit = "ebic")$coefficients[-1]
    met_ebic <- TM_check_metrics(coef_ebic, Beta[1,])
    
    # --- ERIC
    coef_eric <- ic_glmnetboundT(data$x, data$y, crit = "eric")$coefficients[-1]
    met_eric <- TM_check_metrics(coef_eric, Beta[1,])
    
    # --- Store metrics
    results <- list(met_theo, met_cv, met_aic, met_bic, met_ebic, met_eric)
    for (j in 1:length(methods)) {
      acc[j] <- acc[j] + as.numeric(results[[j]]$accuracy) / ntest
      U_fit[j] <- U_fit[j] + as.numeric(results[[j]]$U_fit) / ntest
      O_fit[j] <- O_fit[j] + as.numeric(results[[j]]$O_fit) / ntest
    }
  }
  
  # Output as list of metrics
  return(list(
    Accuracy = 100*acc,
    Underfitted = 100*U_fit,
    Overfitted = 100*O_fit
  ))
}

