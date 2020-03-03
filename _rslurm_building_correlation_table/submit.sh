#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=1
#SBATCH --job-name=building_correlation_table
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
