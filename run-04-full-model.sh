#!/bin/sh
#SBATCH --job-name="04-full-model-right-prior"
#SBATCH --output="04-full-model-right-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.0/lib64/R/bin/Rscript --vanilla \
  simulations/04-full-model/mcmc_right_prior_estimates.r \
	simulations/04-full-model/mcmc_right_prior_estimates.rout &
