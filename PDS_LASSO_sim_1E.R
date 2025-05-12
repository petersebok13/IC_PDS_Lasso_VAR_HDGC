###########################  1_E   #############################################
#1.E Irrepresentable violated, Identity cov matrix
set.seed(42)

n_sims = 1000
sample_sizes = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
sparsity_1E = 0.2
true_driver_1E = sparsity_1E*data_x_dimension
true_driver_1E
sparse_check_1E <- array(0, dim = c(4,4))
dimnames(sparse_check_1E) <- list(data_x_dimension, sample_sizes)
for (i in 1:4) {
  for (j in 1:4) { 
    if (sample_sizes[j] > true_driver_1E[i]*log(data_x_dimension)[j]) {
      sparse_check_1E[i,j] = 1
    }              #sparsity check
  }
}  
sparse_check_1E



beta_min_1E = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_1E[i,j] = sqrt(true_driver_1E[i] * log (data_x_dimension[i]) / sample_sizes[j])
  }
}
beta_min_1E

#simulation

accuracy_1E = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1E)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1E

U_fit_1E = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1E)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1E

O_fit_1E = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1E)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1E

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1E[k], min = 2, max = 2), 
             rep(0, data_x_dimension[k]-true_driver_1E[k]))# selecting betas higher than the required minimum and zero betas
    sigma = 0.5
    Sigma = diag(data_x_dimension[k])
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    if (any(eigen(Sigma)$value <= 0)) {stop("Sigma is not positive definite")}
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1E[k,,j] <- result$Accuracy
      U_fit_1E[k,,j] <- result$Underfitted
      O_fit_1E[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1E);print("Underfit");print(U_fit_1E);print("Overfit"); print(O_fit_1E)

###########################  1_E_2   ###########################################
#1.E Irrepresentable violated, Toeplitz cov matrix

accuracy_1E_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1E_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1E_2

U_fit_1E_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1E_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1E_2

O_fit_1E_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1E_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1E_2

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1E[k], min = 2,max = 2), 
             rep(0, data_x_dimension[k]-true_driver_1E[k])) # selecting betas higher than the required minimum and zero betas
    sigma = 0.5
    Sigma = 0.2^toeplitz(0:(data_x_dimension[k]-1))
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    if (any(eigen(Sigma)$value <= 0)) {stop("Sigma is not positive definite")}
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1E_2[k,,j] <- result$Accuracy
      U_fit_1E_2[k,,j] <- result$Underfitted
      O_fit_1E_2[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}

