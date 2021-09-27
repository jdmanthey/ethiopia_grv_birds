#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cat_vcf
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-39

chr_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

gunzip ${chr_array}.gz

# run vcftools with SNP and invariant site output, 20% max missing data, no indels
vcftools --vcf ${chr_array} --max-missing 0.8 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out ${chr_array%.g.vcf}

bgzip ${chr_array%.g.vcf}.recode.vcf

tabix -p vcf ${chr_array%.g.vcf}.recode.vcf.gz
