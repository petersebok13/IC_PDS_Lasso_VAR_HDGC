###########################  2_E   #############################################
#2.E Irrepresentable cond. violated, Identity cov matrix
set.seed(42)

n_sims = 1000
time_horizons = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_2E = 0.2*data_x_dimension
true_driver_2E
sparse_check_2E <- array(0, dim = c(4,4))
dimnames(sparse_check_2E) <- list(data_x_dimension, time_horizons)
for (i in 1:4) {
  for (j in 1:4) { 
    if (time_horizons[j] > true_driver_2E[i]*log(data_x_dimension)[j]) {
      sparse_check_2E[i,j] = TRUE
    }              #sparsity check
  }
}  
sparse_check_2E

beta_min_2E = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_2E[i,j] = sqrt(true_driver_2E[i] * log (data_x_dimension[i]) / time_horizons[j])
  }
}
beta_min_2E



accuracy_2E = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_2E)= list(data_x_dimension, model_names_basic, n = time_horizons)
accuracy_2E

U_fit_2E = array(data=0, dim = c(4,6,4))
dimnames(U_fit_2E)= list(data_x_dimension, model_names_basic, n = time_horizons)
U_fit_2E

O_fit_2E = array(data=0, dim = c(4,6,4))
dimnames(O_fit_2E)= list(data_x_dimension, model_names_basic, n = time_horizons)
O_fit_2E

for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.45  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:(true_driver_2E[k]+1)) {
        sign_value = ifelse(s %% 2 == 1, 0.45, -0.45)   # odds at +0.45, evens at -0.45
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
      
    }
    sigma = 0.5
    Sigma = diag(data_x_dimension[k]+1)
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #100 simulations
    result <- TM_score_test_VAR(Beta, time_horizons[j], n_sims, lag = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_2E[k,,j] <- result$Accuracy
      U_fit_2E[k,,j] <- result$Underfitted
      O_fit_2E[k,,j] <- result$Overfitted
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_2E);print("Underfit");print(U_fit_2E);print("Overfit"); print(O_fit_2E)

###########################  2_E_2   ###########################################
#2.E Irrepresentable cond. violated, Toeplitz cov matrix

accuracy_2E_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_2E_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
accuracy_2E_2

U_fit_2E_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_2E_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
U_fit_2E_2

O_fit_2E_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_2E_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
O_fit_2E_2

for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.45  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:(true_driver_2E[k]+1)) {
        sign_value = ifelse(s %% 2 == 1, 0.45, -0.45)  # odds at +0.45, evens at -0.45
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
      
    }
    sigma = 0.5
    Sigma = 0.2^toeplitz(0:data_x_dimension[k])
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #100 simulations
    result <- TM_score_test_VAR(Beta, time_horizons[j], n_sims, lag = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_2E_2[k,,j] <- result$Accuracy
      U_fit_2E_2[k,,j] <- result$Underfitted
      O_fit_2E_2[k,,j] <- result$Overfitted
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_2E_2);print("Underfit");print(U_fit_2E_2);print("Overfit"); print(O_fit_2E_2)

