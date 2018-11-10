#!/bin/sh
#SBATCH --job-name="02-missing-wrong-prior"
#SBATCH --output="02-missing-wrong-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.4/lib64/R/bin/Rscript --vanilla \
	simulations/02-missing/mcmc_wrong_prior_estimates.r
