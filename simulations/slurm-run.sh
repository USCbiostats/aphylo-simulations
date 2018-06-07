#!/bin/sh
#SBATCH --job-name=aphylo
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4
#SBATCH --ntasks-per-core=1
SLURM_SUBMIT_DIR = /home/rcf-proj/pdt/vegayon/aphylo-simulations/simulations
#SBATCH --test-only

R CMD BATCH dgp.r dgp.rout