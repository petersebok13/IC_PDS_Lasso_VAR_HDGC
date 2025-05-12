###########################  3A   #############################################
#3.A Conditions hold (Beta_min, irrespresentable, sparsity), Identity cov matrix
set.seed(42)

n_sims = 1000
time_horizons = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_3A = 0.2*data_x_dimension
true_driver_3A
sparse_check_3A <- array(0, dim = c(4,4))
dimnames(sparse_check_3A) <- list(data_x_dimension, time_horizons)
for (i in 1:4) {
  for (j in 1:4) { 
    if (time_horizons[j] > true_driver_3A[i]*log(data_x_dimension)[j]) {
      sparse_check_3A[i,j] = TRUE
    }              #sparsity check
  }
}  
sparse_check_3A

beta_min_3A = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_3A[i,j] = sqrt(true_driver_3A[i] * log (data_x_dimension[i]) / time_horizons[j])
  }
}
beta_min_3A



Power_3A = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(Power_3A)= list(data_x_dimension, model_names_basic, n = time_horizons)
Power_3A

Size_3A = array(data=0, dim = c(4,6,4))
dimnames(Size_3A)= list(data_x_dimension, model_names_basic, n = time_horizons)
Size_3A



for (j in 1:length(time_horizons)) {
  for (k in 1:(length(data_x_dimension))) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.45  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:true_driver_3A[k]) {
        sign_value = ifelse(s %% 2 == 1, 0.45, -0.45)   # odds at +0.45, evens at -0.45
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
      
    }
    
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- HDGC_power_size_score_test(Beta, time_horizons[j], 
                                         n_sims, lag = 1, varcov = "Identity")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Power", "Size") %in% names(result))) {
      Power_3A[k,,j] <- result$Power
      Size_3A[k,,j] <- result$Size

      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Power");print(Power_3A);print("Size");print(Size_3A)

###########################  3_A_2   ###########################################
#3.A Conditions hold (Beta_min, irrespresentable, sparsity), Toeplitz cov matrix

Power_3A_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(Power_3A_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
Power_3A_2

Size_3A_2 = array(data=0, dim = c(4,6,4))
dimnames(Size_3A_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
Size_3A_2


for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.45  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:(true_driver_3A[k])) {
        sign_value = ifelse(s %% 2 == 1, 0.45, -0.45)   # odds at +0.45, evens at -0.45
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
    }
    
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- HDGC_power_size_score_test(Beta, time_horizons[j], 
                                         n_sims, lag = 1, varcov = "Toeplitz")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Power", "Size") %in% names(result))) {
      Power_3A_2[k,,j] <- result$Power
      Size_3A_2[k,,j] <- result$Size

      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}

print("Power");print(Power_3A_2);print("Size");print(Size_3A_2)

print(Power_3A);print(Power_3A_2)
