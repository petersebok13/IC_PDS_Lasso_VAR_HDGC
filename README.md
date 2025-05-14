# IC_PDS_Lasso_VAR_HDGC
This repository accompanies my Master's thesis, which investigates the optimization of the Post-Double Selection (PDS) Lasso-VAR model using various model selection techniques to estimate Granger causality networks in high-dimensional time series settings. 
The repository contains all R scripts necessary to reproduce the simulation studies, methodological comparisons, and empirical analysis conducted in the thesis. The implemented methods include common criteria (AIC, BIC, CV) and advanced approaches (ERIC, EBIC, Theoretical), evaluated under multiple data-generating processes.

Contents (Branches):

simulation/: Scripts to run Monte Carlo simulations for model evaluation

  Part 1: Cross-sectional data and model accuracy,
  Part 2: Time series data and model accuracy,
  Part 3: Time series data and Granger causality

  Scenarios:
  A: All conditions hold,
  B: All conditions fail to hold,
  C: Beta-min condition violated,
  D: Sparsity condition violated,
  E: Irrepresentable condition violated

methods/: All the functions necessary for implementing data generation, model estimation, model selection, score tests (accuracy, power and size), and network estimation

empirical/: Code and the data to replicate the financial volatility network analysis

Sebők, Péter (2025): Information Criterion Based Optimization of the PDS-Lasso-VAR Model and its Application on the High-Dimensional Granger Network Estimation. Corvinus University of Budapest, Institute of Finance.
