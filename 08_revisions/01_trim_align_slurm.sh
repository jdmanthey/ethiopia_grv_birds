#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

module load intel java bwa samtools

# define main working directory
workdir=/lustre/scratch/jmanthey/02_ethiopia_popgen

# define the reference genome
refgenome=/home/jmanthey/references/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna

# run bbduk
/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${SLURM_ARRAY_TASK_ID}_R1.fastq.gz in2=${workdir}/00_fastq/${SLURM_ARRAY_TASK_ID}_R2.fastq.gz out1=${workdir}/01_cleaned/${SLURM_ARRAY_TASK_ID}_R1.fastq.gz out2=${workdir}/01_cleaned/${SLURM_ARRAY_TASK_ID}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bwa mem
bwa mem -t 12 ${refgenome} ${workdir}/01_cleaned/${SLURM_ARRAY_TASK_ID}_R1.fastq.gz ${workdir}/01_cleaned/${SLURM_ARRAY_TASK_ID}_R2.fastq.gz > ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.sam

# convert sam to bam
samtools view -b -S -o ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.bam ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.sam

# remove sam
rm ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.sam

# clean up the bam file
/lustre/work/jmanthey/gatk-4.1.0.0/gatk CleanSam -I ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.bam -O ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}.bam

# sort the cleaned bam file
/lustre/work/jmanthey/gatk-4.1.0.0/gatk SortSam -I ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned.bam -O ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted.bam --SORT_ORDER coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned.bam

# add read groups to sorted and cleaned bam file
/lustre/work/jmanthey/gatk-4.1.0.0/gatk AddOrReplaceReadGroups -I ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted.bam -O ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${SLURM_ARRAY_TASK_ID}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
/lustre/work/jmanthey/gatk-4.1.0.0/gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 -M ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_markdups_metric_file.txt -I ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted_rg.bam -O ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_cleaned_sorted_rg.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${SLURM_ARRAY_TASK_ID}_final.bam
