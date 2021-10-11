#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=gen_phase
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-1170

source activate whatshap

# define reference 
reference=/home/jmanthey/references/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna


# define output directory
outdir=/lustre/scratch/jmanthey/03_ethiopia_popgen/04_msmc


# use individual for this job in the helper array
ind_array=$( head -n${SLURM_ARRAY_TASK_ID} eth_msmc_individual_array.txt | tail -n1 )


# use scaffold for this job in the helper array
scaf_array=$( head -n${SLURM_ARRAY_TASK_ID} eth_msmc_scaffold_array.txt | tail -n1 )


# define some variables
bamfile=${ind_array}
bamfile_full=/lustre/scratch/jmanthey/03_ethiopia_popgen/01_bam_files/${bamfile}
ind=${bamfile%_final.bam}
scaffold=${scaf_array}


# calculate mean coverage for this individual and this scaffold
meancov=`samtools depth -r $scaffold $bamfile_full | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
# echo to log file
echo "Mean coverage for individual $ind, scaffold ${scaffold}: $meancov"


# designate output mask and vcf file names, and phasing stats
mask_ind=${outdir}/ind_mask.${ind}.${scaffold}.samtools.bed.gz

vcf=${outdir}/${ind}.${scaffold}.samtools.vcf

phased_vcf=${outdir}/${ind}.${scaffold}.whatshap.vcf.gz

phased_stats=${outdir}/${ind}.${scaffold}.whatshap.stats.tsv


# run the genotyping and mask-creating script
bcftools mpileup -Ou -r ${scaffold} --threads 4 -f $reference $bamfile_full | bcftools call -c --threads 10 -V indels | ~/msmc-tools-master/bamCaller.py $meancov $mask_ind > ${vcf}


# phase with whatshap and run stats of output
whatshap phase --reference $reference --ignore-read-groups -o $phased_vcf $vcf $bamfile_full

whatshap stats --tsv=$phased_stats $phased_vcf 

