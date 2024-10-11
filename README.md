Supporting Information for
"Developing and Comparing Four Families of Bayesian Network Autocorrelation Models for Binary Outcomes: Estimating Peer Effects Involving Adoption of Medical Technologies"

This file contains 8 scripts for the simulation study of the paper on network autocorrelation model for binary outcomes. 
Each script forms the standalone file with data generation and model estimation for 4 different models (i.e., probit latent response model, logistic latent response model, probit mean model, logistic mean model) using Bayesian approach with different prior distributions of rho.
In particular, scripts ending with "_U" represent estimating the models using uniform prior of rho and scripts ending with "_TU" represent estimating the models using trasnformed uniform prior of rho. 
Robustness analysis is conducted by generating the data under one model and estimating another model (Analysis can be done using the combinations of existing scripts with 16 combinations in total for each prior distribution of rho) and examples of the robustness analysis are under the folder "robustness". 
