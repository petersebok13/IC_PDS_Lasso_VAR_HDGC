#' Simplified 1-lag HDGC with post-double selection and flexible penalty
HDGC_x_y <- function(x, y, data, std = FALSE, alpha = 1, sign = 0.05,
                     method = c("ic", "theoretical", "cv"),
                     crit = c("aic", "bic", "ebic", "eric"),
                     fallback_ridge = TRUE, ridge_lambda = 1e-2) {
  
  method <- match.arg(method)
  crit <- match.arg(crit)
  
  if (!is.data.frame(data)) stop("Data must be a dataframe")
  if (std) data <- as.data.frame(scale(data))
  
  # Step 1: Lagged data (1 lag only)
  lagged_data <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
  for (m in 1:ncol(data)) {
    lagged_data[2:nrow(data), m] <- data[1:(nrow(data) - 1), m]
  }
  colnames(lagged_data) <- paste0(colnames(data), "_lag1")
  
  df <- cbind(data[2:nrow(data), ], lagged_data[2:nrow(data), ])
  rownames(df) <- NULL
  
  # Step 2: Define variables
  y_vec <- df[[y]]
  x_vec <- df[[x]]
  y_lag <- df[[paste0(y, "_lag1")]]
  x_lag <- df[[paste0(x, "_lag1")]]
  Z <- as.matrix(df[, grep("_lag1$", names(df))])
  Z_no_x <- Z[, !colnames(Z) %in% paste0(x, "_lag1"), drop = FALSE]
  
  # Step 3: Post-double selection
  get_selected <- function(X, Y) {
    if (method == "ic") {
      model <- ic_glmnetboundT(X, Y, crit = crit, alpha = alpha, standardize = FALSE)
      return(names(which(model$coefficients != 0)))
    } else if (method == "theoretical") {
      sigma_sq <- mean(lm(Y ~ X)$residuals^2)
      lambda <- sqrt(sigma_sq) / sqrt(nrow(X)) * qnorm(1 - 0.05 / (2 * ncol(X) * log(nrow(X))))
      fit <- glmnet(X, Y, alpha = alpha, lambda = lambda, standardize = FALSE)
      return(rownames(coef(fit))[coef(fit)[, 1] != 0])
    } else if (method == "cv") {
      fit <- rolling_cv_glmnet(X, Y, alpha = alpha, train_ratio = 0.8)
      coef_vec <- coef(fit$final_model)
      selected <- rownames(coef_vec)[coef_vec[, 1] != 0]
      
      return(selected)
    }
  }
  
  s_y <- get_selected(Z_no_x, y_vec)
  s_x <- get_selected(Z_no_x, x_lag)                    
  selected <- as.character(union(s_y, s_x))
  selected <- setdiff(selected, "(Intercept)")
  
  if (length(selected) > 0) {
    S_mat <- Z_no_x[, selected, drop = FALSE]
  } else {
    S_mat <- NULL
  }
  
  # Step 4: Run Granger test model
  test_data <- data.frame(y = y_vec, y_lag = y_lag, x_lag = x_lag)
  if (!is.null(S_mat)) {
    test_data <- cbind(test_data, S_mat)
  }
  
  full_formula <- as.formula(paste("y ~ y_lag + x_lag",
                                   if (!is.null(S_mat)) paste("+", paste(colnames(S_mat), collapse = "+")) else ""))
  restr_formula <- as.formula(paste("y ~ y_lag",
                                    if (!is.null(S_mat)) paste("+", paste(colnames(S_mat), collapse = "+")) else ""))
  
  fit_full <- lm(full_formula, data = test_data)
  fit_restr <- lm(restr_formula, data = test_data)
  
  r2_full <- summary(fit_full)$r.squared
  r2_restr <- summary(fit_restr)$r.squared
  k <- length(coef(fit_full)) - length(coef(fit_restr))  # number of restrictions (x_lag)
  df1 <- k
  df2 <- nrow(test_data) - length(coef(fit_full))
  
  fstat <- ((r2_full - r2_restr) / df1) / ((1 - r2_full) / df2)
  pval <- pf(fstat, df1, df2, lower.tail = FALSE)
  
  if (is.na(pval)) {
    warning(paste0("NA p-value for test: ", x, " at ", y))
    return(1)
  } else if (pval <= sign) {
    cat(paste0(x, " Granger-causes ", y, ", p-value = ", signif(pval, 4), "\n"))
    return(pval)
  } else {
    cat(paste0(x, " does NOT Granger-cause ", y, ", p-value = ", signif(pval, 4), "\n"))
    return(pval)
  }
}
