#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=25
#SBATCH --job-name=sponge_interaction_filter_singlesubset
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
