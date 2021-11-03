#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=bed_mask
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

cd /lustre/scratch/jmanthey/03_ethiopia_popgen/04_genotype_mask

ind_array=$( head -n${SLURM_ARRAY_TASK_ID} individual_list.txt | tail -n1 )

ind_array=${ind_array}.genotyped.sites

while read -r scaffold1; do
        gunzip -cd ${SLURM_ARRAY_TASK_ID}.coverage.mask | grep $scaffold1 > ${ind_array%.genotyped.sites}.${scaffold1}.coverage.mask
        
        gzip ${ind_array%.genotyped.sites}.${scaffold1}.coverage.mask
        
done < scaffold_list.txt

