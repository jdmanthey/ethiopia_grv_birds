#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cat_vcf
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-39

chr_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

bgzip ${chr_array}

tabix -p vcf ${chr_array}.gz
