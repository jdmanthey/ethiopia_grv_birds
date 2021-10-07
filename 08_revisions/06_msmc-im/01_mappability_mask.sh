
conda create --name py27 python=2.7
source activate py27
conda install -c anaconda zlib
module load gnu7/7.3.0
# install the snpable package

# interactively run the snpable pipeline
interactive -p quanah -N 1 -c 24

# activate python 2 environment
source activate py27

# navigate to working directory
cd /lustre/scratch/jmanthey
mkdir mappability
cd mappability

# set up reference and kmer number variables
reference=/home/jmanthey/references/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna
k=35

# run splitfa and split into multiple files
/home/jmanthey/seqbility-20091110/splitfa $reference $k | split -l 20000000

# cat the output together
cat x* >> ref_split.${k}

# run bwa aln to align kmers to reference
bwa aln -t 24 -R 100000 -O 3 -E 3 ${reference} ref_split.${k} > ref_split.${k}.sai

# create sam output from sai file
bwa samse -f ref_split.${k}.sam $reference ref_split.${k}.sai ref_split.${k}

# generate raw mask from sam file
/home/jmanthey/seqbility-20091110/gen_raw_mask.pl ref_split.${k}.sam > ref_split.${k}.prelim_mask.fa

# generate final mask file from raw mask
/home/jmanthey/seqbility-20091110/gen_mask ref_split.${k}.prelim_mask.fa > ref_split.${k}.final_mask.fa


