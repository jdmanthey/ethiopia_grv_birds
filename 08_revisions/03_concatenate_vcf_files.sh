#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cat_vcf
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-25

chr_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_cat.txt | tail -n1 )

chr_array2=${chr_array%__}

chr_start=${chr_array}__a.g.vcf

chr_output=${chr_array}.g.vcf

grep "#" $chr_start > $chr_output

for i in $( ls ${chr_array}__*vcf ); do grep -v "#" $i >> $chr_output; done
