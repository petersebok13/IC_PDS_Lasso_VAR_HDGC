###########################  1_A   #############################################
#1.A Conditions hold (Beta_min, irrespresentable, sparsity), Identity cov matrix
set.seed(42)

n_sims = 1000
sample_sizes = c(50,100,200,500)

#sparsity ( n > s*log(p) )
data_x_dimension = c(10,20,50,100)
true_driver_1A = 0.2*data_x_dimension
true_driver_1A
sparse_check_1A <- array(0, dim = c(4,4))
dimnames(sparse_check_1A) <- list(data_x_dimension, sample_sizes)
for (i in 1:4) {
 for (j in 1:4) { 
  if (sample_sizes[j] > true_driver_1A[i]*log(data_x_dimension)[j]) {
    sparse_check_1A[i,j] = TRUE
  }              #sparsity check
 }
}  
sparse_check_1A

beta_min_1A = array(data=0, dim = c(4,4))
for (i in 1:4) {
  for (j in 1:4) {
    beta_min_1A[i,j] = sqrt(true_driver_1A[i] * log (data_x_dimension[i]) / sample_sizes[j])
  }
}
beta_min_1A



accuracy_1A = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1A)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1A

U_fit_1A = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1A)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1A

O_fit_1A = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1A)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1A

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1A[k], min = 2, max = 2), 
             rep(0, data_x_dimension[k]-true_driver_1A[k])) # selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Identity")

    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1A[k,,j] <- result$Accuracy
      U_fit_1A[k,,j] <- result$Underfitted
      O_fit_1A[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1A);print("Underfit");print(U_fit_1A);print("Overfit"); print(O_fit_1A)

###########################  1_A_2   ###########################################
#1.A Conditions hold (Beta_min, irrespresentable, sparsity), Toeplitz cov matrix

accuracy_1A_2 = array(data=0, dim = c(4,6,4))
model_names_basic <- c("Theo","CV", "AIC", "BIC", "EBIC", "ERIC")
dimnames(accuracy_1A_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
accuracy_1A_2

U_fit_1A_2 = array(data=0, dim = c(4,6,4))
dimnames(U_fit_1A_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
U_fit_1A_2

O_fit_1A_2 = array(data=0, dim = c(4,6,4))
dimnames(O_fit_1A_2)= list(data_x_dimension, model_names_basic, n = sample_sizes)
O_fit_1A_2

for (j in 1:length(sample_sizes)) {
  for (k in 1:length(data_x_dimension)) {
    
    #sampling of betas meeting the criteria
    beta = c(runif(true_driver_1A[k], min = 2,max = 2), 
             rep(0, data_x_dimension[k]-true_driver_1A[k])) # selecting betas higher than the required minimum and zero betas
    
    #100 simulations
    result <- TM_score_test_basic_v2(beta, sample_sizes[j], n_sims, sigma = 1, varcov = "Toeplitz")
    
    # Ensure `result` is in the correct format (named numeric vector or list)
    if (is.list(result) && 
        all(c("Accuracy", "Underfitted", "Overfitted") %in% names(result))) {
      accuracy_1A_2[k,,j] <- result$Accuracy
      U_fit_1A_2[k,,j] <- result$Underfitted
      O_fit_1A_2[k,,j] <- result$Overfitted
      cat("Sample size:", sample_sizes[j], "- p:", data_x_dimension[k], "\n")
    } else {
      stop("Unexpected output format.")
    }
  }
}
print("Accuracy");print(accuracy_1A_2);print("Underfit");print(U_fit_1A_2);print("Overfit"); print(O_fit_1A_2)


