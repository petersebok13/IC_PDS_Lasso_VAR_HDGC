# IC_PDS_Lasso_VAR_HDGC
This repository accompanies my Master's thesis, which investigates the optimization of the Post-Double Selection (PDS) Lasso-VAR model using various model selection techniques to estimate Granger causality networks in high-dimensional time series settings. 
The repository contains all R scripts necessary to reproduce the simulation studies, methodological comparisons, and empirical analysis conducted in the thesis. The implemented methods include common criteria (AIC, BIC, CV) and advanced approaches (ERIC, EBIC, Theoretical), evaluated under multiple data-generating processes.

Contents:

simulation/: Scripts to run Monte Carlo simulations for model evaluation

methods/: Functions implementing different model selection techniques

empirical/: Code to replicate the financial volatility network analysis

utils/: Helper functions used across simulations and applications

plots/: Reproducible figures used in the thesis

Sebők, Péter (2025): Information Criterion Based Optimization of the PDS-Lasso-VAR Model and its Application on the High-Dimensional Granger Network Estimation
