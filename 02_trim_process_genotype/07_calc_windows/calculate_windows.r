# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# file name for parsing chromosome and window bounds
filename_simple <- strsplit(args[1], "windows/")[[1]][2]

# add functions for calculations
source("window_stat_calculations.r")

# define minimum number of sites to keep a fasta file 
min_sites <- 50000

# no scientific notation
options(scipen=999)

# read in input file
input_file <- read.table(args[1], sep="\t", stringsAsFactors=F)

# subset input file 
input_file_genotypes <- input_file[,4:ncol(input_file)]

# read in populations
populations <- read.table(args[2], sep="\t", stringsAsFactors=F, header=T)

# define output name
output_name <- paste(strsplit(args[1], ".simple")[[1]][1], "__stats.txt", sep="")

# write output file
write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=output_name, sep="\t")

# calculate heterozygosity for each individual
heterozygosity(input_file, populations, output_name, filename_simple)

# calculate differentiation statistics for each pairwise comparison
all_combinations <- combn(unique(populations$PopName), 2)
# keep only those where the species is the same in the combination
all_combinations <- all_combinations[,sapply(strsplit(all_combinations[1,], "_"), "[[", 1) == sapply(strsplit(all_combinations[2,], "_"), "[[", 1)]
for(a in 1:ncol(all_combinations)) {
	# define populations
	a_pop1 <- all_combinations[1,a]
	a_pop2 <- all_combinations[2,a]
	

	# subset vcf inputs
	a_input1 <- input_file_genotypes[,populations$PopName == a_pop1]
	if(is.null(dim(a_input1))) { a_input1 <- as.matrix(a_input1) }
	# remove sites that are not either invariant or bi-allelic SNPs
	a_input1 <- a_input1[nchar(input_file[,2]) == 1 & nchar(input_file[,3]) == 1, ]
	if(is.null(dim(a_input1))) { a_input1 <- as.matrix(a_input1) }

	a_input2 <- input_file_genotypes[,populations$PopName == a_pop2]
	if(is.null(dim(a_input2))) { a_input2 <- as.matrix(a_input2) }
	# remove sites that are not either invariant or bi-allelic SNPs
	a_input2 <- a_input2[nchar(input_file[,2]) == 1 & nchar(input_file[,3]) == 1, ]
	if(is.null(dim(a_input2))) { a_input2 <- as.matrix(a_input2) }

	# determine which sites have missing data and remove them
	a_input <- cbind(a_input1, a_input2)
	dont_keep <- apply(a_input, 1, any_missing)
	a_input1 <- a_input1[dont_keep == FALSE, ]
	if(is.null(dim(a_input1))) { a_input1 <- as.matrix(a_input1) }
	a_input2 <- a_input2[dont_keep == FALSE, ]
	if(is.null(dim(a_input2))) { a_input2 <- as.matrix(a_input2) }
			
	differentiation(a_input1, a_input2, a_pop1, a_pop2, output_name, filename_simple)
	
}
