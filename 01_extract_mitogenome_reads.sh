#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N eth_six_mt
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-30

/lustre/work/jmanthey/bbmap/bbsplit.sh \
in1=/lustre/scratch/jmanthey/02_ethiopia_popgen/00_fastq/${SGE_TASK_ID}_R1.fastq.gz \
in2=/lustre/scratch/jmanthey/02_ethiopia_popgen/00_fastq/${SGE_TASK_ID}_R2.fastq.gz \
ref=/lustre/scratch/jmanthey/02_ethiopia_popgen/bird_mitogenomes.fasta \
basename=/lustre/scratch/jmanthey/02_ethiopia_popgen/01_no_mtdna/${SGE_TASK_ID}%.fastq.gz \
outu1=/lustre/scratch/jmanthey/02_ethiopia_popgen/01_no_mtdna/${SGE_TASK_ID}_R1.fastq.gz \
outu2=/lustre/scratch/jmanthey/02_ethiopia_popgen/01_no_mtdna/${SGE_TASK_ID}_R2.fastq.gz
