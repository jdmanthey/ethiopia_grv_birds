# use a simplified vcf as input for creating a structure file with snps at least minimum_dist apart

# column headers
vcf_colnames <- c("POS", "REF", "ALT", "1", "2", "3", "4", "5")
options(scipen=999)

	# set up names of individuals
	ind_names <- c("Cossypha_semirufa_EB008", "Cossypha_semirufa_EB009", "Cossypha_semirufa_EB024", "Cossypha_semirufa_EB044", "Cossypha_semirufa_EB062")


# minimum distance between snps
minimum_dist <- 20000


# order individuals west to east
individuals_ordered <- ind_names[c(4,5,1,2,3)]

# all the files to read
x_files <- list.files(pattern="*struc.simple.vcf")

# determine unique chromosomes
x_unique <- unique(sapply(strsplit(x_files, "__"), "[[", 1))

# loop for each chromosome
for(a in 1:length(x_unique)) {
	print(a)
	print(paste("Reading files"))
	# read in all the files and concatenate them
	b_files <- x_files[sapply(strsplit(x_files, "__"), "[[", 1) %in% x_unique[a]]
	for(b in 1:length(b_files)) {
		if(b == 1) {
			a_rep <- read.table(b_files[b], stringsAsFactors=F, sep="")
		} else {
			a_rep <- rbind(a_rep, read.table(b_files[b], stringsAsFactors=F, sep=""))
		}
	}
	print(paste("Removing indels"))
	# remove indels
	a_rep <- a_rep[nchar(a_rep[,2]) == 1 & nchar(a_rep[,3]) == 1, ]
	print(paste("Removing variants with missing data"))
	# remove sites with missing data
	for(b in 4:ncol(a_rep)) {
		a_rep <- a_rep[grepl(pattern="\\./\\.", a_rep[,b]) == F,]
	}
	print(paste("Subsetting SNPs by distance"))
	# site positions and distances
	site_positions <- a_rep[,1]
	keep <- site_positions[1]
	site_positions <- site_positions[site_positions > site_positions[1] + minimum_dist]
	while(length(site_positions) > 0) {
		keep <- c(keep, site_positions[1])
		site_positions <- site_positions[site_positions > site_positions[1] + minimum_dist]
	}
	new_a_rep <- a_rep[match(keep, a_rep[,1]),]
	
	# modify all genotypes to indicate the alleles rather than 0/1
	for(b in 4:ncol(new_a_rep)) {
		# remove phasing information
		new_a_rep[,b] <- gsub("\\|", "/", new_a_rep[,b])
		
		# replace 0 and 1 with arbitrary letters
		new_a_rep[,b] <- gsub("0", "w", new_a_rep[,b])
		new_a_rep[,b] <- gsub("1", "v", new_a_rep[,b])		
	}
	for(b in 1:nrow(new_a_rep)) {
		# replace those arbitrary letters with the correct genotypes for structure
		new_a_rep[b, 4:ncol(new_a_rep)] <- gsub("w", new_a_rep[b,2], new_a_rep[b, 4:ncol(new_a_rep)])
		new_a_rep[b, 4:ncol(new_a_rep)] <- gsub("v", new_a_rep[b,3], new_a_rep[b, 4:ncol(new_a_rep)])		
	}
	for(b in 4:ncol(new_a_rep)) {
		# replace genotype letters with numbers for structure
		new_a_rep[,b] <- gsub("A", "1", new_a_rep[,b])
		new_a_rep[,b] <- gsub("C", "2", new_a_rep[,b])
		new_a_rep[,b] <- gsub("G", "3", new_a_rep[,b])
		new_a_rep[,b] <- gsub("T", "4", new_a_rep[,b])
	}
	
	
	# write names to column headers
	name_header <- c()
	for(b in 1:length(ind_names)) {
		name_header <- c(name_header, ind_names[b])
		name_header <- c(name_header, ind_names[b])
	}
	if(a == 1) {
		write(name_header, file="cossypha_20000.structure", sep="\t", ncolumns=50)
	}

	# prepare the matrix for writing to file
	output <- list()
	for(b in 4:ncol(new_a_rep)) {
		output[[((b-3)*2 - 1)]] <- as.character(sapply(strsplit(new_a_rep[,b], "/"), "[[", 1))
		output[[((b-3)*2)]] <- as.character(sapply(strsplit(new_a_rep[,b], "/"), "[[", 2))
	}
	# combine the list
	output2 <- c()
	for(b in 1:length(output)) {
		if(b == 1) {
			output2 <- output[[1]]
		} else {
			output2 <- cbind(output2, output[[b]])
		}
	}
	# write to output
	write.table(output2, file="cossypha_20000.structure", sep="\t", quote=F, row.names=F, col.names=F, append=T)

}

# rewrite the output transposed
x <- read.table(file="cossypha_20000.structure", sep="\t", stringsAsFactors=F)
x <- t(x)
write.table(x, file="cossypha_20000b.structure", sep="\t", quote=F, row.names=F, col.names=F)


