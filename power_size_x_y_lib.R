#########
#Functions

HDGC_power_size_score_test <- function(Beta, n, ntest, lag = 1, varcov = "Identity") {
  
  methods <- c("Theo", "CV", "AIC", "BIC", "EBIC", "ERIC")
  Power <- setNames(rep(0, length(methods)), methods)
  Size  <- setNames(rep(0, length(methods)), methods)
  
  for (i in 1:ntest) {
    
    data <- gen_sim_data_ts_GC(Beta_matrix = Beta, n = n, varcov = varcov, lag = lag)
    data_df <- data.frame(data)
    colnames(data_df) <- paste0("V", 1:ncol(data))
    
    # Rename second column and last column
    colnames(data_df)[ncol(data_df)] <- "Size_Scope"   
    colnames(data_df)[2] <- "Power_Scope"
    
    for (method in methods) {
      
      if (method == "Theo") {
        crit <- NULL
        penalty <- "theoretical"
      } else if (method == "CV") {
        crit <- NULL
        penalty <- "cv"
      } else {
        crit <- tolower(method)   # "aic", "bic", "ebic", "eric"
        penalty <- "ic"
      }
      
      ## Size Scope Test
      result_size <- HDGC_x_y(x = "Size_Scope", y = "V1", data = data_df, std = FALSE,
                              #lags = lag, 
                              crit = crit, alpha = 1, sign = 0.05,
                              method = penalty)
      
      ## Power Scope Test
      result_power <- HDGC_x_y(x = "Power_Scope", y = "V1", data = data_df, std = FALSE,
                               #lags = lag, 
                               crit = crit, alpha = 1, sign = 0.05,
                               method = penalty)
      
      # Update counters
      if (result_size < 0.05) {
        Size[method] <- Size[method] + 1
      }
      
      if (result_power < 0.05) {
        Power[method] <- Power[method] + 1
      }
    }
  }
  
  # Average over the number of simulations
  Power <- 100*Power / ntest
  Size  <- 100*Size / ntest
  
  return(list(
    Power = Power,
    Size = Size
  ))
}
