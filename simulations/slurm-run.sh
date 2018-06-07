#!/bin/bash
#SBATCH --job-name=aphylo
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4
#SBATCH --ntasks-per-core=1
##SLURM_SUBMIT_DIR = /home/rcf-proj/pdt/vegayon/aphylo-simulations/simulations
##SBATCH --test-only

# Relevant environment variables
export RPROJECT=/home/rcf-proj/pdt/vegayon/aphylo-simulations
export PANTHER=/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1\_altVersion/hmmscoring/PANTHER13.1
export R_LIBS=~/R/x86_64-redhat-linux-gnu-library/3.4/

/usr/usc/R/3.4.0/bin/Rscript --vanilla dgp.r
