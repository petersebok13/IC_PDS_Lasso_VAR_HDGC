###########################  3B   #############################################
#3.B Conditions do NOT hold (Beta_min, irrespresentable, sparsity), Identity cov matrix
set.seed(42)

n_sims = 1000
time_horizons = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_3B = 0.8*data_x_dimension
true_driver_3B
sparse_check_3B <- array(0, dim = c(4,4))
dimnames(sparse_check_3B) <- list(data_x_dimension, time_horizons)
for (i in 1:4) {
  for (j in 1:4) { 
    if (time_horizons[j] > true_driver_3B[i]*log(data_x_dimension)[j]) {
      sparse_check_3B[i,j] = TRUE
    }              #sparsity check
  }
}  
sparse_check_3B

beta_min_3B = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_3B[i,j] = sqrt(true_driver_3B[i] * log (data_x_dimension[i]) / time_horizons[j])
  }
}
beta_min_3B



Power_3B = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(Power_3B)= list(data_x_dimension, model_names_basic, n = time_horizons)
Power_3B

Size_3B = array(data=0, dim = c(4,6,4))
dimnames(Size_3B)= list(data_x_dimension, model_names_basic, n = time_horizons)
Size_3B



for (j in 1:length(time_horizons)) {
  for (k in 1:(length(data_x_dimension))) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.15  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:true_driver_3B[k]) {
        sign_value = ifelse(s %% 2 == 1, 0.15, -0.15)  # odds at +0.15, evens at -0.15
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
      
    }
    sigma = 0.5
    Sigma = diag(data_x_dimension[k]+1)
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- HDGC_power_size_score_test(Beta, time_horizons[j], 
                                         n_sims, lag = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Power", "Size") %in% names(result))) {
      Power_3B[k,,j] <- result$Power
      Size_3B[k,,j] <- result$Size
      
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Power");print(Power_3B);print("Size");print(Size_3B)

###########################  3_B_2   ###########################################
#3.B Conditions do NOT hold (Beta_min, irrespresentable, sparsity), Toeplitz cov matrix

Power_3B_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(Power_3B_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
Power_3B_2

Size_3B_2 = array(data=0, dim = c(4,6,4))
dimnames(Size_3B_2)= list(data_x_dimension, model_names_basic, n = time_horizons)
Size_3B_2


for (j in 1:length(time_horizons)) {
  for (k in 1:length(data_x_dimension)) {
    
    Beta = diag(data_x_dimension[k] + 1) * 0.15  # +1 because of y
    
    for (i in 1:nrow(Beta)) {
      
      Beta_row = rep(0, ncol(Beta))  # initialize the full row with zeros
      
      # Insert alternating signs (+, -, +, - ...)
      for (s in 1:(true_driver_3B[k])) {
        sign_value = ifelse(s %% 2 == 1, 0.15, -0.15)  # odds at +0.15, evens at -0.15
        Beta_row[s] = sign_value
      }
      
      # Assign the new row into Beta
      Beta[i,] = Beta_row
    }
    sigma = 0.5
    Sigma = 0.2^toeplitz(0:data_x_dimension[k])
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    #if (any(abs(eigen(Beta)$value)) >= 1) { stop("Beta is not compatible.") }
    
    #1000 simulations
    result <- HDGC_power_size_score_test(Beta, time_horizons[j], 
                                         n_sims, lag = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Power", "Size") %in% names(result))) {
      Power_3B_2[k,,j] <- result$Power
      Size_3B_2[k,,j] <- result$Size
      
      cat("Sample size:", time_horizons[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Power");print(Power_3B_2);print("Size");print(Size_3B_2)

print(Power_3B);print(Power_3B_2)
