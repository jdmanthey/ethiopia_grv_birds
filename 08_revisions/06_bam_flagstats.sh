#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

module load intel java bwa samtools

# define main working directory
workdir=/lustre/scratch/jmanthey/03_ethiopia_popgen

# get mapping stats for each final bam file
samtools flagstat ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_final.bam > ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_flagstat.txt


