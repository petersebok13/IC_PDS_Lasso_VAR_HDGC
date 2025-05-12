###########################  2_C   #############################################
#2.C Beta-min violated, Identity cov matrix
set.seed(42)

n_sims = 1000
time_horizons = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_2C = 0.2*data_x_dimension
true_driver_2C
sparse_check_2C <- array(0, dim = c(4,4))
dimnames(sparse_check_2C) <- list(data_x_dimension, time_horizons)
for (i in 1:4) {
  for (j in 1:4) { 
    if (time_horizons[j] > true_driver_2C[i]*log(data_x_dimension)[j]) {
      sparse_check_2C[i,j] = TRUE
    }              #sparsity check
  }
}  
sparse_check_2C

beta_min_2C = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_2C[i,j] = sqrt(true_driver_2C[i] * log (data_x_dimension[i]) / time_horizons[j])
  }
}
beta_min_2C



accuracy_2C = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_2C)= list(data_x_dimension, model_names_basic, n = time_horizons)
accuracy_2C

U_fit_2C = array(data=0, dim = c(4,6,4))
dimnames(U_fit_2C)= list(data_x_dimension, model_names_basic, n = time_horizons)
U_fit_2C

O_fit_2C = array(data=0, dim = c(4,6,4))
dimnames(O_fit_2C)= list(data_x_dimension, model_names_basic, n = time_horizons)
O_fit_2C

for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.095  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:true_driver_2C[k]) {
        sign_value = ifelse(s %% 2 == 1, 0.095, -0.095)  # odds at +0.095, evens at -0.095
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
      
    }
    
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- TM_score_test_VAR(Beta, time_horizons[j], n_sims, lag = 1, varcov = "Identity")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_2C[k,,j] <- result$Accuracy
      U_fit_2C[k,,j] <- result$Underfitted
      O_fit_2C[k,,j] <- result$Overfitted
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_2C);print("Underfit");print(U_fit_2C);print("Overfit"); print(O_fit_2C)

###########################  2_C_2   ###########################################
#2.C Beta-min violated, Toeplitz cov matrix

accuracy_2C_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_2C_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
accuracy_2C_2

U_fit_2C_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_2C_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
U_fit_2C_2

O_fit_2C_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_2C_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
O_fit_2C_2

for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.095  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:(true_driver_2C[k])) {
        sign_value = ifelse(s %% 2 == 1, 0.095, -0.095)  # odds at +0.095, evens at -0.095
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
    }
    
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- TM_score_test_VAR(Beta, time_horizons[j], n_sims, lag = 1, varcov = "Toeplitz")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_2C_2[k,,j] <- result$Accuracy
      U_fit_2C_2[k,,j] <- result$Underfitted
      O_fit_2C_2[k,,j] <- result$Overfitted
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_2C_2);print("Underfit");print(U_fit_2C_2);print("Overfit"); print(O_fit_2C_2)

