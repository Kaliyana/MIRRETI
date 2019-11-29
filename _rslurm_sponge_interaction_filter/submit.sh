#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=40
#SBATCH --job-name=sponge_interaction_filter
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
