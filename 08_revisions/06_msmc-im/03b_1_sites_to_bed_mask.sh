#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

work_dir=/lustre/scratch/jmanthey/03_ethiopia_popgen/03_vcf/

cd $work_dir

ind_array=$( head -n${SLURM_ARRAY_TASK_ID} individual_list.txt | tail -n1 )

out_dir=/lustre/scratch/jmanthey/03_ethiopia_popgen/04_genotype_mask/

# filter the vcfs for one individual each with no missing data
for i in $( ls *g.vcf ); do
        
        vcftools --vcf $i --indv ${ind_array} --max-missing 1 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out ${out_dir}${i%.g.vcf}_${ind_array}

done

# change to the fitlered vcf directory
cd $out_dir

# cat all the site information to a new file
for i in $( ls *_${ind_array}.recode.vcf ); do

        grep -v "#" $i | cut -f 1,2 >> ${ind_array}.genotyped.sites
        
done

