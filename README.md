Supporting Information for
"Developing and Comparing Four Families of Bayesian Network Autocorrelation Models for Binary Outcomes: Estimating Peer Effects Involving Adoption of Medical Technologies"

This file contains R scripts for the simulation study of the paper "Developing and Comparing Four Families of Bayesian Network Autocorrelation Models for Binary Outcomes: Estimating Peer Effects Involving Adoption of Medical Technologies". 
Each script forms the standalone file with data generation and model estimation for 4 different models (i.e., probit latent response model, logistic latent response model, probit mean model, logistic mean model) using Bayesian approach with different prior distributions of rho.
In particular, scripts ending with "_U" represent estimating the models using uniform prior of rho and scripts ending with "_TU" represent estimating the models using trasnformed uniform prior of rho. 
Robustness analysis is conducted by generating the data under one model and estimating another model (Analysis can be done using the combinations of existing scripts with 16 combinations in total for each prior distribution of rho) and examples of the robustness analysis are under the folder "robustness". For illustration, the true value of rho is set to 0 in each script, and these scripts could be easily adapted to other true values of rho (i.e., -0.2 and 0.2).  

robust_plot.R is the R script for generating the robust plot: Figure 1 in the paper using the simulation results saved in robust_dat.xlsx. 

The real data used for the motivating analyses contains patient identifiable information and so cannot be made available. However, the network size, density and distribution of edge weights used in the simulation are very similar to the network in the motivating example and therefore can be easily adapted to analyze a real data set. We also generate a synthetic dataset (SynthDS.RData) for users to analyze to emulate the real data analyses using the probit latent response model as an example: SynthDS_analysis.R.
