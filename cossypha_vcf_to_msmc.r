# use a simplified vcf as input for creating msmc files (needs to include variant and invariant sites)

# create msmc files for each individual, one per chromosome
library(data.table)

# define some variables
out_dir <- "cossypha_demography"

# column headers
vcf_colnames <- c("POS", "REF", "ALT", "1", "2", "3", "4", "5")
options(scipen=999)

# set up names of individuals
ind_names <- c("C_semirufa_EB008", "C_semirufa_EB009", "C_semirufa_EB024", "C_semirufa_EB044", "C_semirufa_EB062")

# all the files to read
x_files <- list.files(pattern="*msmc.simple.vcf")

# functions to be applied with apply across rows for each allele
# col1 = position, col2 = reference allele, col3 = alt alleles, col4 = genotype
replace_genotype <- function(xxx) {
	if(xxx[4] == "0/0") { 
		alleles <- c(xxx[2], xxx[2])
	} else {
		potential_alleles <- c(xxx[2], strsplit(xxx[3], ",")[[1]])
		alleles <- c(potential_alleles[as.numeric(strsplit(xxx[4], "/")[[1]][1]) + 1],
						potential_alleles[as.numeric(strsplit(xxx[4], "/")[[1]][2]) + 1])
	}
	return(alleles)
}


# determine unique chromosomes
x_unique <- unique(sapply(strsplit(x_files, "__"), "[[", 1))

# make output directory
dir.create(out_dir)

# loop for each chromosome
for(a in 1:length(x_unique)) {
	writeLines(paste("Reading files"))
	# read in all the files and concatenate them
	b_files <- x_files[sapply(strsplit(x_files, "__"), "[[", 1) %in% x_unique[a]]
	for(b in 1:length(b_files)) {
		writeLines(paste("Reading file ", b, " for scaffold ", a, sep=""))
		if(b == 1) {
			a_rep <- read.table(b_files[b], stringsAsFactors=F, sep="")
		} else {
			a_rep <- rbind(a_rep, read.table(b_files[b], stringsAsFactors=F, sep=""))
		}
	}
	# loop for each individual for this chromosome
	for(b in 1:length(ind_names)) {
		writeLines(paste("Individual", b))
		# make output directory for this individual
		if(a == 1) { dir.create(paste(out_dir, "/", ind_names[b], sep="")) }
		# subset locations and variants for this individual
		b_rep <- a_rep[,c(1,2,3,b+3)]
		b_rep <- b_rep[b_rep[,4] != "./.",]
		# replace all phased genotypes to keep it simple
		b_rep[,4] <- gsub("\\|", "/", b_rep[,4])
		
		# change row names of sites kept so far
		rownames(b_rep) <- seq(from=1, to=nrow(b_rep), by=1)
		
		# subset only variants
		b_rep_matches <- b_rep[,4] == "0/0" | b_rep[,4] == "1/1" | b_rep[,4] == "2/2" | b_rep[,4] == "3/3" | b_rep[,4] == "4/4" | 
							b_rep[,4] == "5/5" | b_rep[,4] == "6/6" 
		b_rep_matches <- b_rep_matches == FALSE
		b_rep_variant <- b_rep[b_rep_matches, ]
		
		# replace genotype with actual variants to remove null alleles from calculations
		b_alleles <- t(as.matrix(apply(b_rep_variant, 1, replace_genotype)))
		b_alleles <- cbind(b_rep_variant[,1], b_alleles)
		b_alleles2 <- data.frame(site=as.numeric(b_alleles[,1]), allele1=as.character(b_alleles[,2]), allele2=as.character(b_alleles[,3]), total_sites=as.numeric(rownames(b_alleles)))
		b_alleles2 <- b_alleles2[grepl("\\*", b_alleles2[,2]) == FALSE, ]
		b_alleles2 <- b_alleles2[grepl("\\*", b_alleles2[,3]) == FALSE, ]
		
		# determine number of sites before each match
		b_rep_sites <- c(0, b_alleles2$total_sites)
		b_rep_sites <- diff(b_rep_sites)
		
		# output
		output <- data.frame(chr=as.character(paste("scaf", a, sep="")), site=as.numeric(b_alleles2[,1]), g_sites=as.numeric(b_rep_sites), genotype=as.character(paste(b_alleles2[,2], b_alleles2[,3], sep="")))
		
		# write output
		write.table(output, file=paste(out_dir, "/", ind_names[b], "/", output[1,1], ".txt", sep=""), sep="\t",
							quote=F, row.names=F, col.names=F)
	}
}




