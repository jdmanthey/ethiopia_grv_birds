#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=gzip
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

gzip ${SLURM_ARRAY_TASK_ID}.coverage.mask
