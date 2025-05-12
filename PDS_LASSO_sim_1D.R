###########################  1_D   #############################################
#1.D Beta min violated, Identity cov matrix
set.seed(42)

n_sims = 1000
sample_sizes = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_1D = 0.2*data_x_dimension
true_driver_1D
sparse_check_1D <- array(0, dim = c(4,4))
dimnames(sparse_check_1D) <- list(data_x_dimension, sample_sizes)
for (i in 1:4) {
  for (j in 1:4) { 
    if (sample_sizes[j] > true_driver_1D[i]*log(data_x_dimension)[j]) {
      sparse_check_1D[i,j] = TRUE
    }              #sparsity check
  }
}  
sparse_check_1D

beta_min_1D = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_1D[i,j] = sqrt(true_driver_1D[i] * log (data_x_dimension[i]) / sample_sizes[j])
  }
}
beta_min_1D



accuracy_1D = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1D)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1D

U_fit_1D = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1D)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1D

O_fit_1D = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1D)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1D

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1D[k], min = 0.095, max = 0.095), 
             rep(0, data_x_dimension[k]-true_driver_1D[k])) # selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Identity")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1D[k,,j] <- result$Accuracy
      U_fit_1D[k,,j] <- result$Underfitted
      O_fit_1D[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1D);print("Underfit");print(U_fit_1D);print("Overfit"); print(O_fit_1D)

###########################  1_D_2   ###########################################
#1.D Beta min violated, Toeplitz cov matrix

accuracy_1D_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1D_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1D_2

U_fit_1D_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1D_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1D_2

O_fit_1D_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1D_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1D_2

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1D[k], min = 0.095, max = 0.095), 
             rep(0, data_x_dimension[k]-true_driver_1D[k])) # selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Toeplitz")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1D_2[k,,j] <- result$Accuracy
      U_fit_1D_2[k,,j] <- result$Underfitted
      O_fit_1D_2[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1D_2);print("Underfit");print(U_fit_1D_2);print("Overfit"); print(O_fit_1D_2)
