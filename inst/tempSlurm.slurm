#!/bin/sh
#SBATCH --job-name=WalkingDead2_simulations
#SBATCH --time=12:00:00
#SBATCH --account=kingspeak
#SBATCH --partition=phs7045

module load R/4.2.2 imagemagick/7.0.8-7
Rscript --vanilla tempSlurm.R
