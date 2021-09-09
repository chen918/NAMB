This file contains 6 scripts for the simulation study of the paper on binary network autocorrelation model. 
Each script forms the standalone file with data generation and model estimation using Bayesian approach with different prior distributions of rho.
In particular, BNAM_uniform_1 uses a uniform prior of rho with a range from 1/lamda_min to 1/lamda_max; BNAM_uniform_2 uses a uniform prior of rho with a range from -1 to 1; BNAM_uniform_3 uses a uniform prior of rho with a range from rho_true-0.25 to rho_true+0.25; BNAM_uniform_T uses a transformed uniform prior of rho.
BNAM_JR uses a Jeffrey rule prior of rho; BNAM_IJ uses an Independence Jeffreys prior.
For each script, the output includes the posterior median estimator of coffefients and rho. In adiition, it also includes the mean squared error of rho, 95% coverage rate of rho, 95% equal-tailed credible interval of rho and average width of credible intrerval of rho. 
