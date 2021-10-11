options(scipen=999)

# read in genome index and only keep scaffolds at least 1 Mbp
genome_index <- read.table("GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna.fai", stringsAsFactors=F)
genome_index <- genome_index[genome_index[,2] >=1000000,]

# remove sex chromosome
genome_index2 <- genome_index[genome_index[,1] != "NC_011493.1",]
write(genome_index2[,1], file="scaffold_list.txt", ncolumns=1)

# read in list of all bam files with format {individualname}_final.bam
bam_file_list <- read.table("ethiopia_bam_list.txt", stringsAsFactors=F)

# loop to make helper files 
individual_array <- c()
scaffold_array <- c()
for(a in 1:nrow(genome_index)) {
	a_scaffold <- genome_index[a,1]
	for(b in 1:nrow(bam_file_list)) {
		b_individual <- bam_file_list[b,1]
		individual_array <- c(individual_array, b_individual)
		scaffold_array <- c(scaffold_array, a_scaffold)
	}
}

write(individual_array, file="eth_msmc_individual_array.txt", ncolumns=1)
write(scaffold_array, file="eth_msmc_scaffold_array.txt", ncolumns=1)
