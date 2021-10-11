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

# haplotype order for each species for east
hap_order_east=("0,1,2,3,4,5" "0,1,2,3" "0,1,2,3,4,5" "0,1,2,3" "0,1,2,3,4,5" "0,1,2,3")

# haplotype order for each species for west
hap_order_west=("6,7,8,9" "4,5,6,7,8,9" "6,7,8,9,10,11" "4,5" "6,7,8,9,10,11" "4,5,6,7,8,9")

# choose the variables for this array
species=${species_array[$array_slot]}
het=${het_array[$array_slot]}
hap_east=${hap_order_east[$array_slot]}
hap_west=${hap_order_west[$array_slot]}


# run the program for the east
/home/jmanthey/msmc2_linux64bit -I ${hap_east} -i 20 -t 18 -m ${het} -s -p 1*2+20*1+1*2+1*3 -o ${species}/east ${species}*multihetsep.txt

# run the program for the west
/home/jmanthey/msmc2_linux64bit -I ${hap_west} -i 20 -t 18 -m ${het} -s -p 1*2+20*1+1*2+1*3 -o ${species}/west ${species}*multihetsep.txt

