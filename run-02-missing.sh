#!/bin/sh
#SBATCH --job-name="02-missing-right-prior"
#SBATCH --output="02-missing-right-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.0/lib64/R/bin/Rscript --vanilla \
	simulations/02-missing/mcmc_right_prior_estimates.r \
	simulations/02-missing/mcmc_right_prior_estimates.rout 
