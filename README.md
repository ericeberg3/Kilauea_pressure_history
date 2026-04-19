Kilauea Pressure History Modeling
This repository contains MATLAB scripts for estimating Kilauea volcano's pressure history and volume changes using Markov Chain Monte Carlo (MCMC) inversions of geodetic data.
HOW TO RUN THE INVERSIONS:
Step 1: Generate Prior Probability Distributions. Before running any MCMC inversions, you must generate the prior probability distribution data. The main MCMC algorithms depend on this data to evaluate the prior probability of proposed model parameters.
	•	Ensure that the required posterior distribution images are located in the Data/post_im/ directory (e.g., x_hmm_kyle.png, vol_hmm.png).
	•	Open MATLAB and run the script: prior_prob.m
	•	This script reads the images, fits appropriate statistical distributions (Normal, Lognormal, Gamma, Uniform) to the parameters via kernel density estimation, and saves the output to a file named Data/paramDists.mat.
Step 2: Prepare Geodetic Data Ensure your GPS and InSAR data vectors, along with their respective covariance matrices, are prepared.
	•	Use the scripts in the Make_cov_insar/ directory (Prep_InSAR_cov.m) to build the data covariance matrices for the ascending and descending InSAR datasets.
	•	Use create_MCMC_data.m (or create_MCMC_data_vol.m) to assemble your final data structures.
Step 3: Configure the MCMC Sampler (mcmc.m), you need to set up your initial states:
	•	Define your initial parameter estimates (x0 or minit), step sizes (xstep), and bounds (xbnds).
Step 4: Execute the MCMC Inversion Run your chosen inversion script, passing in your forward model function, data, initial parameters, and the required number of iterations.
	•	Run mcmc.m. This will use the Metropolis-Hastings algorithm and output the accepted samples, log-likelihoods, and acceptance ratios.
Step 5: Analyze and Visualize Results Once the chain has finished running:
	•	Discard the initial "burn-in" period of the generated Markov chains.
	•	Navigate to the Plotting functions/ directory.
	•	Use scripts like plotHists.m to view parameter posteriors, plot_insar.m or plot_insar_new.m to compare observed versus modeled surface displacements, and makeplots.m for general diagnostic plotting.
