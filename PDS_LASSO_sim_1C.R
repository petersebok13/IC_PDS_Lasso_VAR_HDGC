###########################  1_C   #############################################
#1.C sparsity violated, Identity cov matrix
set.seed(42)

n_sims = 1000
sample_sizes = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
sparsity_1C = 0.8
true_driver_1C = sparsity_1C*data_x_dimension
true_driver_1C
sparse_check_1C <- array(0, dim = c(4,4))
dimnames(sparse_check_1C) <- list(data_x_dimension, sample_sizes)
for (i in 1:4) {
  for (j in 1:4) { 
    if (sample_sizes[j] > true_driver_1C[i]*log(data_x_dimension)[j]) {
      sparse_check_1C[i,j] = 1
    }              #sparsity check
  }
}  
sparse_check_1C



beta_min_1C = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_1C[i,j] = sqrt(true_driver_1C[i] * log (data_x_dimension[i]) / sample_sizes[j])
  }
}
dimnames(beta_min_1C) <- list(d = data_x_dimension,n = sample_sizes)
beta_min_1C

#simulation

accuracy_1C = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1C)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1C

U_fit_1C = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1C)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1C

O_fit_1C = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1C)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1C

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1C[k], min = 3, max = 3), 
             rep(0, data_x_dimension[k]-true_driver_1C[k]))# selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Identity")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1C[k,,j] <- result$Accuracy
      U_fit_1C[k,,j] <- result$Underfitted
      O_fit_1C[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1C);print("Underfit");print(U_fit_1C);print("Overfit"); print(O_fit_1C)

###########################  1_C_2   ###########################################
#1.C Sparsity violated, Toeplitz cov matrix

accuracy_1C_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1C_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1C_2

U_fit_1C_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1C_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1C_2

O_fit_1C_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1C_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1C_2

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1C[k], min = 3,max = 3), 
             rep(0, data_x_dimension[k]-true_driver_1C[k])) # selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Toeplitz")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1C_2[k,,j] <- result$Accuracy
      U_fit_1C_2[k,,j] <- result$Underfitted
      O_fit_1C_2[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1C_2);print("Underfit");print(U_fit_1C_2);print("Overfit"); print(O_fit_1C_2)


