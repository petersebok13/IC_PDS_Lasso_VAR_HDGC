#################################
#basic functions
gen_sim_data <- function(beta, sample_size, varcov = "Identity", sigma = 1) {
  
  p <- length(beta)  # Number of predictors
  
  if (is.character(varcov) && varcov == "Identity") {
    
    x <- rnorm_multi(
      n = sample_size, 
      mu = rep(0, p),
      sd = rep(1, p),
      r = rep(0, p * (p - 1) / 2),  # Identity matrix at 0 correlations
      empirical = FALSE, 
      as.matrix = TRUE
    )
    
  } else if (is.character(varcov) && varcov == "Toeplitz") {
    
    # Toeplitz correlation matrix
    varcov_matrix <- 0.5^toeplitz(0:(p - 1))
    
    # Convert to lower triangular vector for `rnorm_multi()`
    varcov_vector <- varcov_matrix[lower.tri(varcov_matrix)]
    
    x <- rnorm_multi(
      n = sample_size, 
      mu = rep(0, p),
      sd = rep(1, p),
      r = varcov_vector,  # Use extracted lower-triangle values
      empirical = FALSE, 
      as.matrix = TRUE
    )
  } else { #check this option: intended for the irrepresentable condition
    
    varcov_vector <- varcov[lower.tri(varcov)]
    
    x <- rnorm_multi(
      n = sample_size, 
      mu = rep(0, p),
      sd = rep(1, p),
      r = varcov_vector,  # Use extracted lower-triangle values
      empirical = FALSE, 
      as.matrix = TRUE
    )
  }
  
  # Generate error term and response variable
  eps <- rnorm(n = sample_size, mean = 0, sd = sigma)
  y <- x %*% beta + eps  # True Data-Generating Process (DGP)
  
  return(list(x = x, y = y))
}

################################

gen_sim_data_ts <- function(Beta_matrix, n, lag = 1, varcov = "Identity") {
  library(tsDyn)
  
  p <- nrow(Beta_matrix) # number of variables
  
  # Build covariance matrix
  if (is.character(varcov)) {
    if (varcov == "Identity") {
      Sigma <- diag(p)
    } else if (varcov == "Toeplitz") {
      Sigma <- 0.5^toeplitz(0:(p-1))
    } else {
      stop("varcov must be 'Identity' or 'Toeplitz'")
    }
  } else {
    # User-specified full variance-covariance matrix
    Sigma <- varcov
  }
  
  # Simulate VAR
  sim <- tsDyn::VAR.sim(B = Beta_matrix, n = n, include = "none", varcov = Sigma)
  Y <- as.matrix(sim)
  
  # Build lagged design matrix
  lagged_data <- embed(Y, lag + 1)
  X <- lagged_data[, -(1:p)]      # remove current time t columns
  Y_target <- lagged_data[, 1]     # model first variable, can generalize later
  
  return(list(x = X, y = Y_target))
}
###############

gen_sim_data_ts_GC <- function(Beta_matrix, n, lag = 1, varcov = "Identity") {
  library(tsDyn)
  
  p <- nrow(Beta_matrix) # number of variables
  
  # Build covariance matrix
  if (is.character(varcov)) {
    if (varcov == "Identity") {
      Sigma <- diag(p)
    } else if (varcov == "Toeplitz") {
      Sigma <- 0.5^toeplitz(0:(p-1))
    } else {
      stop("varcov must be 'Identity' or 'Toeplitz'")
    }
  } else {
    # User-specified full variance-covariance matrix
    Sigma <- varcov
  }
  
  # Simulate VAR
  sim <- tsDyn::VAR.sim(B = Beta_matrix, n = n, include = "none", varcov = Sigma)
  Y <- as.data.frame(sim)
  

  return(Data = Y)
}

###############
rolling_cv_glmnet <- function(x, y, alpha = 1, lambda_grid = NULL, train_ratio = 0.8) {
  library(glmnet)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # Default lambda grid if not provided
  if (is.null(lambda_grid)) {
    lambda_grid <- exp(seq(log(0.001), log(10), length.out = 250))  #CHECK
  }
  
  n_train <- floor(train_ratio * n)
  n_valid <- n - n_train
  
  # Split the data
  x_train <- x[1:n_train, , drop = FALSE]
  y_train <- y[1:n_train]
  
  x_valid <- x[(n_train + 1):n, , drop = FALSE]
  y_valid <- y[(n_train + 1):n]
  
  # Fit the model on training data for all lambdas
  fit <- glmnet(x_train, y_train, alpha = alpha, lambda = lambda_grid)
  
  # Predict on validation set
  preds <- predict(fit, newx = x_valid)
  
  # Calculate Mean Squared Error for each lambda
  mse <- colMeans((preds - y_valid)^2)
  
  # Find best lambda (smallest validation error)
  best_idx <- which.min(mse)
  best_lambda <- lambda_grid[best_idx]
  
  # Refit final model on full training set with best lambda
  final_model <- glmnet(x_train, y_train, alpha = alpha, lambda = best_lambda)
  
  return(list(
    final_model = final_model,
    best_lambda = best_lambda,
    validation_mse = mse,
    lambda_grid = lambda_grid
  ))
}



