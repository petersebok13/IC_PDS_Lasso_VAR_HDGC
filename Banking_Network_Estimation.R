#### Empirical network estimation of financial institutions ####
library(readxl)
#setwd
logreturn = read_excel("financial_data.xlsx", sheet = "LogReturns")
volatility = logreturn[,-1]^2
#View(logreturn)
#View(volatility)

#Networks

#Theoretical
NetGC_HDGC(volatility, log = FALSE, std = TRUE,
           alpha = 1, sign = 0.05, crit = "aic",
           lags = 1, method = "theoretical",
           plot = "veins", cluster = FALSE, verbose = TRUE,
           title = "Theoretical rule")

#CV
NetGC_HDGC(volatility, log = FALSE, std = TRUE,
           alpha = 1, sign = 0.05, crit = "aic",
           lags = 1, method = "cv",
           plot =  "veins", cluster = FALSE, verbose = TRUE,
           title = "CV rule")

#AIC
NetGC_HDGC(volatility, log = FALSE, std = T,
                       alpha = 1, sign = 0.05, crit = "aic",
                       lags = 1, method = "ic",
                       plot = "veins", cluster = FALSE, verbose = TRUE,
           title = "AIC rule")

#BIC
NetGC_HDGC(volatility, log = FALSE, std = T,
           alpha = 1, sign = 0.05, crit = "bic",
           lags = 1, method = "ic",
           plot = "veins", cluster = FALSE, verbose = TRUE,
           title = "BIC rule")

#EBIC
NetGC_HDGC(volatility, log = FALSE, std = T,
           alpha = 1, sign = 0.05, crit = "ebic",
           lags = 1, method = "ic",
           plot =  "veins", cluster = FALSE, verbose = TRUE,
           title = "EBIC rule")

#ERIC
NetGC_HDGC(volatility, log = FALSE, std = T,
           alpha = 1, sign = 0.05, crit = "eric",
           lags = 1, method = "ic",
           plot =  "veins", cluster = FALSE, verbose = TRUE,
           title = "ERIC rule")
