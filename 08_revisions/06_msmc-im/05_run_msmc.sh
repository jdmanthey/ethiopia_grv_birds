#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=msmc_all
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=18
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-6

# array slot
array_slot=$(($SLURM_ARRAY_TASK_ID - 1))

# define species
species_array=(cossypha crithagra melaenornis parophasma turdus zosterops)

# define half the heterozygosity (mean per genus) for the MSMC r parameter
het_array=(0.00229627 0.00132413 0.00139479 0.00043941 0.00215914 0.00147503)

# define allele array since not all samples have the same sample sizes
allele_array=(0,0,0,0,1,1,1,1 0,0,0,0,1,1,1,1 0,0,0,0,1,1,1,1 0,0,0,0,1,1 0,0,0,0,1,1,1,1 0,0,0,0,1,1,1,1)

# choose the variables for this array
species=${species_array[$array_slot]}
het=${het_array[$array_slot]}
allele=${allele_array[$array_slot]}

mkdir ${species}

# run the program for the east
/home/jmanthey/msmc_1.1.0_linux64bit -P ${allele} -i 20 -t 18 --fixedRecombination --skipAmbiguous -p 1*2+20*1+1*2+1*3 -o ${species}/msmc ${species}*multihetsep.txt

