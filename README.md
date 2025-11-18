# MOO

This repository contains an R implementation of the MOO (masking-one-out) algorithms. The main idea of the MOO procedure is to mask an observed entry and compare it with its imputed values. We implement proposed MOO, MOORT (masking-one-out rank transformation), and MOOEN (masking-one-out with energy distance) in the simulation and data application. See more details in our paper: https://arxiv.org/abs/2511.10048.

The Simulation folder includes the functions for implementing the proposed algorithms, functions for the imputers (mean imputation, Expectation-Maximization (EM) algorithm, nearest-neighbor hot deck (NN HD), complete-case missing value (CCMV), Markov missing graph (MMG), and multiple imputation by chained equations (MICE), and helper functions for G-MMG and EM. The EM, CCMV, and MMG methods fit Gaussian models, where MMG uses the Gaussian-MMG specification. We provide examples in MOO.R and combined.R to illustrate how to implement the algorithms and how to visualize PI diagram. 

The Read Data folder provides additional functions for implementing variable-wise MOO, MOORT, and MOOEN, as well as functions for the imputers and helper functions for MMG.  
