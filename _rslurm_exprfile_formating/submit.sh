#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --job-name=exprfile_formating
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
