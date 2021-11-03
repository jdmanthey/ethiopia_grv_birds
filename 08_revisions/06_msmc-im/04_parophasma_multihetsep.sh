#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=msmc_zost
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-38

# define individual array for input names (can be numbers or characters, depending on naming scheme)
ind_array=(26 30 27 28)

# read in scaffolds to use
chr=$( head -n${SLURM_ARRAY_TASK_ID} scaffold_list.txt | tail -n1 )

# prefix of output
prefix="zosterops"

# define working directory
work_dir=/lustre/scratch/jmanthey/03_ethiopia_popgen/04_msmc/

# define output directory (make sure it is already made)
out_dir=/lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/

# list of vcf files for each individual
vcf_list=( "${ind_array[@]/%/.${chr}.whatshap.vcf.gz}" )
vcf_list=( "${vcf_list[@]/#/${work_dir}}" )
vcf_list=$(printf " %s" "${vcf_list[@]}")
echo $vcf_list

# list of mask files for each individual
mask_list=( "${ind_array[@]/%/.${chr}.coverage.mask.gz}" )
mask_list=( "${mask_list[@]/#/--mask ${work_dir}}" )
mask_list=$(printf " %s" "${mask_list[@]}")
echo $mask_list

/home/jmanthey/msmc-tools-master/generate_multihetsep.py $mask_list $vcf_list > ${out_dir}${prefix}.${chr}.multihetsep.txt
