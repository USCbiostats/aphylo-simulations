#!/bin/sh
#SBATCH --job-name="04-full-model-wrong-prior"
#SBATCH --output="04-full-model-wrong-prioir-slurm%A.out"
#SBATCH --time="5:00:00"
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
/usr/usc/R/3.4.4/lib64/R/bin/Rscript --vanilla \
  simulations/04-full-model/mcmc_wrong_prior_estimates.r 
