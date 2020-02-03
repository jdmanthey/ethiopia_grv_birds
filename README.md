# ethiopia_grv_birds
Differentiation in 6 species of birds across the Great Rift Valley (GRV)

Methods:

    1. 00_prep_fastq: run the rename.sh script to rename raw fastq files
    2. 01_mitogenomes: use bbsplit to separate putative mitochondrial reads for later use in Geneious to assemble mitogenomes
    3. 02_trim_process_genotype: filter and genotype all samples
        a. 01_filter_process.sh: filter the sequences and create bam alignment files
        b. 02_calc_depth.sh: use samtools to calculate seq. coverage (summarize with 02b script)
        c. 03_create_genotype_scripts.r: R script to make cluster submission scripts for GATK genotyping
        d. 04_create_file_list.sh: list all vcf output for summarization
        e. 05_filter_vcfs.r: R script to make cluster submission scripts for vcf filtering
    4. 03_structure: use each of the R scripts to convert filtered vcf files to STRUCTURE input files
    5. 04_mutation_rate: calculate mutation rates from vcf files
        a. 00_cat_vcf.sh: combine all vcf files for a single vcf
        b. 01_subset_gff.r: subset the Zebra Finch GFF file to only take CDS
        c. 02_bedtools_intersect.sh: use bedtools to extract the CDS regions from the VCF
        d. 03_extract...: R script to filter and write alignment for each gene from the extracted VCF
        e. 04_4d_sites.r: R script to further filter alignments and make 4-fold degenerate site alignments
        f. 05_trim_phylo.r: R script to trim the Passeriformes reference phylogeny to the needed taxa here
        g. 06_phyml...: information and calculatations about running model testing, phylogeny estimate, etc.
    6. 05_demography: demographic anaylses using MSMC
        a. run all the ____vcf_to_msmc.r files to get MSMC input files from the vcf inputs
        b. 02_create_bootstraps.sh: create bootstrap reps with the MSMC utility scripts
        c. 03_create_array_job.r: R script to write cluster submission array job with helper files
        d. 04_plot_demography.r: R script to plot the MSMC output
        e. 05_harmonic_mean_calc.r: R script to calculate harmonic mean pop. sizes from MSMC output
    7. 06_niche_modeling: scripts to run niche modeling
        a. 01_rarefy_points: R script to rarefy the downloaded points from GBIF
        b. 02_modeling_plus_stats: R script to perform ENM and some stats on the output models
