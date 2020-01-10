# function to test if there is any data in a row (use with apply)
test_if_data <- function(xxxx) {
	xxxx <- length(xxxx[xxxx != "./."])
	return(xxxx)
}

# function to return one non-missing base per site (to use with apply for species-specific vcf subset)
allele_return <- function(xxxx) {
	if("0/0" %in% xxxx) {
		return(0)
	} else if("1/1" %in% xxxx) {
		return(1)
	} else {
		xxxx <- xxxx[xxxx != "./."]
		xxxx <- gsub("\\|", "/", xxxx)
		xxxx <- strsplit(xxxx, "/")[[1]][1]
		return(xxxx)
	}
}

# function to replace genotype info with the actual bases for each species
# xxxx = species genotypes list, xxx_list = allele_list
bp_input <- function(xxxx, xxx_list) {
	return(xxx_list[xxxx])
}


# read in gff
gff <- read.table("finch_cds2.gff", sep="\t", stringsAsFactors=F)

# read in vcf
vcf <- read.table("genes_filtered.simple.vcf", sep="\t", stringsAsFactors=F)

# out_directory
dir.create("fasta")

# subset the gff to those chromosomes that were genotyped
gff <- gff[gff[,1] %in% unique(vcf[,1]), ]

# header
vcf_header <- c("CHROM","POS","REF","ALT","1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","23","24","25","26","27","28","29","3","30","4","5","6","7","8","9")
colnames(vcf) <- vcf_header

# read in popmap
popmap <- read.table("eth_popmap.txt", sep="\t", stringsAsFactors=F)
# add species names to popmap
popmap <- cbind(popmap, sapply(strsplit(popmap[,1], "_EB"), "[[", 1))
colnames(popmap) <- c("Individual", "Number", "Species")

# unique cds in gff file
cds <- unique(gff[,9])

# remove highly overlapping cds (same start sites as a conservative take)
cds_test1 <- gff[match(cds, gff[,9]),]
cds_test2 <- paste(cds_test1[,1], cds_test1[,4])
cds <- cds_test1[match(unique(cds_test2), cds_test2), 9]

# filter vcf for sites found in at least one individual per species
keep_list <- list()
for(a in 1:length(unique(popmap[,3]))) {
	writeLines(paste("Species", a))
	# subset vcf for this species
	a_rep <- vcf[, colnames(vcf)[colnames(vcf) %in% popmap[popmap[,3] == unique(popmap[,3])[a], 2]]]
	keep_list[[a]] <- apply(a_rep, 1, test_if_data)
}
# determine which snps are found at least once in each species
keep <- keep_list[[1]] > 0 & keep_list[[2]] > 0 & keep_list[[3]] > 0 & keep_list[[4]] > 0 & keep_list[[5]] > 0 & keep_list[[6]] > 0

# subset the vcf
vcf <- vcf[keep,]

# loop for each unique cds to determine whether to keep
keep <- c()
for(a in 1:length(cds)) {
	if(a %% 50 == 0) {print(a)}
	# subset gff
	a_gff <- gff[gff[,9] == cds[a],]
	
	# keep going if the CDS is divisible by 3 (codon length and full)
	if(sum(a_gff[,5] - a_gff[,4] + 1) %% 3 == 0) {
		# filter vcf for this cds
		a_vcf <- c()
		for(b in 1:nrow(a_gff)) {
			a_vcf <- rbind(a_vcf, vcf[vcf[,1] == a_gff[b,1] & vcf[,2] >= a_gff[b,4] & vcf[,2] <= a_gff[b,5],])
		}
		
		# check if the nrow of this subset is at least 95% the total length of the cds
		if(nrow(a_vcf) >= (0.95 * sum(a_gff[,5] - a_gff[,4] + 1))) {
			keep <- c(keep, a)
		}
	}
}

# garbage collection
gc()

# only analyze those passing above filtering
cds_new <- cds[keep]
# loop for each cds to extract sequence
for(a in 1:length(cds_new)) {
	if(a %% 50 == 0) {print(a)}
	# subset gff
	a_gff <- gff[gff[,9] == cds_new[a],]
	
	# filter vcf for this cds
	a_vcf <- c()
	for(b in 1:nrow(a_gff)) {
		a_vcf <- rbind(a_vcf, vcf[vcf[,1] == a_gff[b,1] & vcf[,2] >= a_gff[b,4] & vcf[,2] <= a_gff[b,5],])
	}
	
	# determine which sites are missing and fill them in with N
	cds_positions <- c()
	for(b in 1:nrow(a_gff)) {
		cds_positions <- c(cds_positions, seq(from=a_gff[b,4], to=a_gff[b,5], by=1))
	}
	cds_missing <- cds_positions[cds_positions %in% a_vcf[,2] == FALSE]
	if(length(cds_missing) > 0) {
		for(b in 1:length(cds_missing)) {
			a_vcf <- rbind(a_vcf, c(a_vcf[1,1], cds_missing[b], "N", ".", rep("0/0", nrow(popmap))))
		}
	}
	
	# order the matrix in the correct order based on orientation of the cds
	if(a_gff[1,7] == "+") {
		a_vcf <- a_vcf[order(a_vcf[,2], decreasing=F),]
	} else {
		a_vcf <- a_vcf[order(a_vcf[,2], decreasing=T),]
	}
	
	# determine alleles for each site
	allele_list <- paste(a_vcf[,3], a_vcf[,4], sep=",")
	# complement the alleles if the cds is in reverse orientation
	if(a_gff[1,7] == "-") {
		allele_list <- gsub("A", 1, allele_list)
		allele_list <- gsub("C", 2, allele_list)
		allele_list <- gsub("G", 3, allele_list)
		allele_list <- gsub("T", 4, allele_list)
		allele_list <- gsub(1, "T", allele_list)
		allele_list <- gsub(2, "G", allele_list)
		allele_list <- gsub(3, "C", allele_list)
		allele_list <- gsub(4, "A", allele_list)
	}
	# replace null alleles with Ns
	allele_list <- gsub("\\*", "N", allele_list)
	allele_list <- strsplit(allele_list, ",")
	
	# alleles for each species to use in alignment
	allele_species <- list()
	for(b in 1:length(unique(popmap[,3]))) {
		# subset vcf for this species
		b_rep <- a_vcf[, colnames(a_vcf)[colnames(a_vcf) %in% popmap[popmap[,3] == unique(popmap[,3])[b], 2]]]
		allele_species[[b]] <- as.numeric(apply(b_rep, 1, allele_return))
		allele_species[[b]] <- paste(mapply(bp_input, as.list(allele_species[[b]] + 1), allele_list), collapse="")
	}
	# add the zebra finch genotype
	b <- b+1
	allele_species[[b]] <- paste(mapply(bp_input, as.list(rep(1, length(allele_list))), allele_list), collapse="")

	
	# write output
	sp_names <- c(as.character(unique(popmap[,3])), "Taeniopygia_guttata")
	sp_names <- paste(">", sp_names, sep="")
	out_name <- paste("fasta/", sapply(strsplit(cds_new[a], "="), "[[", 2), ".fasta", sep="")
	for(b in 1:length(sp_names)) {
		if(b == 1) {
			write(sp_names[b], file=out_name, ncolumns=1)
		} else {
			write(sp_names[b], file=out_name, ncolumns=1, append=T)
		}
			write(allele_species[[b]], file=out_name, ncolumns=1, append=T)
	}
	
}




















