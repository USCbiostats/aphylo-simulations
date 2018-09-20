#!/bin/sh
#SBATCH --job-name="01-gold-standard-right-prior"
#SBATCH --output="01-golds-standard-right-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.0/lib64/R/bin/Rscript --vanilla \
	simulations/01-gold-standard/mcmc_right_prior_estimates.r \
	simulations/01-gold-standard/mcmc_right_prior_estimates.rout
