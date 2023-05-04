#!/bin/sh
#SBATCH --job-name=WalkingDead2_simulations
#SBATCH --time=8:00:00
#SBATCH --account= owner-guest
#SBATCH --partition= notchpeak-guest

module load R/4.2.2
Rscript --vanilla tempSlurm.R
