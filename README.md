# Wolbachia-invasion-dynamics-non-spatial
gyp_sim.cpp (version 04/10/2016) was developed by Dr Penelope A. Hancock

gyp_sim.cpp comprises C++ code for implementing the model
of mosquito-Wolbachia dynamics developed in Hancock et al. 2016,
"Predicting Wolbachia invasion dynamics in Aedes aegypti populations using
models of density-dependent demographic traits", BMC Biology.

Also available on Figshare: doi: 10.6084/m9.figshare.3980472

To run gyp_sim use the command gyp_sim.exe<gyp_sim_inits.txt

The file gyp_sim_inits.txt contains the following inputs:

Cohort_means.txt (a file for storing the mean development times of
uninfected larvae in each cohort)

Cohort_stds.txt (a file for storing the standard deviations of the 
development times of uninfected larvae in each cohort)

mu_p.txt (a file for storing the number of uninfected pupae that eclose
on each day)

L_file.txt (a file for storing the number of uninfected larvae present
on each day)

A_file.txt (a file for storing the number of uninfected adults present
on each day)

lambda.txt (a file for storing the per-capita female fecundity at the time
that each cohort is hatched)

0.07 (additional density-INdependent daily mortality experienced by adults 
in the field environment)

0.01 (additional density-INdependent daily mortality experienced by larvae 
 in the field environment)
