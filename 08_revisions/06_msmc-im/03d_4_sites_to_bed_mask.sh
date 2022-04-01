#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=bed_mask
#SBATCH --nodes=1 --ntasks=8
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-29

module load intel R

cd /lustre/scratch/jmanthey/03_ethiopia_popgen/04_genotype_mask

ind_array=$( head -n${SLURM_ARRAY_TASK_ID} individual_list.txt | tail -n1 )

ind_array=${ind_array}.genotyped.sites

Rscript genotyped_sites_bed.r $ind_array
