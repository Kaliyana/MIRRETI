#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=25
#SBATCH --job-name=runtime_analysis_sponge_6
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
