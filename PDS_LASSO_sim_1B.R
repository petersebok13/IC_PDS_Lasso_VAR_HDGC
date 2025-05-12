###########################  1_B   #############################################
#1.B Conditions not hold (Beta_min, irrespresentable, sparsity), Identity cov matrix
set.seed(42)

n_sims = 1000
sample_sizes = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
sparsity_1B = 0.8
true_driver_1B = sparsity_1B*data_x_dimension
true_driver_1B
sparse_check_1B <- array(0, dim = c(4,4))
dimnames(sparse_check_1B) <- list(data_x_dimension, sample_sizes)
for (i in 1:4) {
  for (j in 1:4) { 
    if (sample_sizes[j] > true_driver_1B[i]*log(data_x_dimension)[j]) {
      sparse_check_1B[i,j] = 1
    }              #sparsity check
  }
}  
sparse_check_1B



beta_min_1B = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_1B[i,j] = sqrt(true_driver_1B[i] * log (data_x_dimension[i]) / sample_sizes[j])
  }
}
beta_min_1B

#simulation

accuracy_1B = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1B)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1B

U_fit_1B = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1B)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1B

O_fit_1B = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1B)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1B

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1B[k], min = 0.15, max = 0.15), 
             rep(0, data_x_dimension[k]-true_driver_1B[k]))# selecting betas higher than the required minimum and zero betas
    sigma = 0.5
    Sigma = diag(data_x_dimension[k])
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    if (any(eigen(Sigma)$value <= 0)) {stop("Sigma is not positive definite")}
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1B[k,,j] <- result$Accuracy
      U_fit_1B[k,,j] <- result$Underfitted
      O_fit_1B[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1B);print("Underfit");print(U_fit_1B);print("Overfit"); print(O_fit_1B)

###########################  1_B_2   ###########################################
#1.B Conditions do not hold (Beta_min, irrespresentable, sparsity), Toeplitz cov matrix

accuracy_1B_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1B_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1B_2

U_fit_1B_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1B_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1B_2

O_fit_1B_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1B_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1B_2

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1B[k], min = 0.15,max = 0.15), 
             rep(0, data_x_dimension[k]-true_driver_1B[k])) # selecting betas higher than the required minimum and zero betas
    sigma = 0.5
    Sigma = 0.2^toeplitz(0:(data_x_dimension[k]-1))
    Sigma = generate_block_cov(Sigma = Sigma, rho = sigma)
    if (any(eigen(Sigma)$value <= 0)) {stop("Sigma is not positive definite")}
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = Sigma)
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1B_2[k,,j] <- result$Accuracy
      U_fit_1B_2[k,,j] <- result$Underfitted
      O_fit_1B_2[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1B_2);print("Underfit");print(U_fit_1B_2);print("Overfit"); print(O_fit_1B_2)


