#!/bin/sh
#SBATCH --job-name="mcmc_right_prior"
#SBATCH --output="mcmc_right_prior%A.out"
#SBATCH --time="02:00:00"
source /usr/usc/R/3.4.0/setup.sh
Rscript --vanilla mcmc_right_prior_estimates.r