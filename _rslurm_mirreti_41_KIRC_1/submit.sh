#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=10
#SBATCH --job-name=mirreti_41_KIRC_1
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
