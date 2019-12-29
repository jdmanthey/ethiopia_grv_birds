#functions used to calculate differentiation between two populations or diversity statistics within populations
# make sure the output file is already written in format (header):
# population1	population2	statistic	chromosome	start	end	number_sites	number_variable_sites	calculated_stat
# e.g., write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=outname, sep="\t")

# input is simplified vcf (entire vcf) and 3-column popmap (ordered the same as the vcf) and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
heterozygosity <- function(xxx, popmap, outname, filename) {
	options(scipen=999)
	xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
	for(a in 1:nrow(popmap)) {
		output_rep <- c(popmap[a,1], "none", "heterozygosity", strsplit(filename, ":")[[1]][1], 
			as.numeric(strsplit(strsplit(filename, ":")[[1]][2], "-")[[1]][1]),
			as.numeric(strsplit(strsplit(filename, "-")[[1]][2], ".simple")[[1]][1]))
		# select this individual
		a_rep <- xxx[,a+3]
		# remove missing genotypes
		a_rep <- a_rep[a_rep != "./."]
		# count number of sites
		a_total <- length(a_rep)
		# remove phasing information
		a_rep <- gsub("\\|", "/", a_rep)
		# find number of heterozygous sites
		a_het <- length(a_rep[a_rep == "0/1"])
		# add to output
		output_rep <- c(output_rep, a_total, a_het, a_het/a_total)
		write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
	}
}


# function to determine if a site is only missing data for a population (used across rows of vcf)
# used in differentiation function
total_missing <- function(xxx) {
	return(length(xxx[xxx != "./."]) > 0)
}	
# function to determine if a site has any missing data (used across rows of vcf)
any_missing <- function(xxx) {
	return(length(xxx[xxx == "./."]) > 0)
}		
# function to determine if a site is polymorphic (to be applied across rows) after removing missing
# used in differentiation function
polymorphic_function2 <- function(xxx) {
	xxx <- xxx[xxx != "./."]
	return(length(unique(xxx)))
}			

# input is two simplified vcfs (subsampled to single population), already filtered for invariant/biallelic SNPs, 
# the names of the populations, and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
# fst is calculation of Reich et al. 2009
# fst is the reich et al. 2009 estimator for small sample sizes
# equation presented nicer in Willing et al. 2012 page 9


differentiation <- function(xxx1, xxx2, popname1, popname2, outname, filename) {
	# remove phasing information
	for(a in 1:ncol(xxx1)) {
		xxx1[,a] <- gsub("\\|", "/", xxx1[,a])
	}
	for(a in 1:ncol(xxx2)) {
		xxx2[,a] <- gsub("\\|", "/", xxx2[,a])
	}
	
	# remove sites that are completely missing from either population
	keep1 <- apply(xxx1, 1, total_missing)
	keep2 <- apply(xxx2, 1, total_missing)
	xxx1 <- xxx1[keep1 == TRUE & keep2 == TRUE, ]
	if(is.null(dim(xxx1))) { xxx1 <- as.matrix(xxx1) }
	xxx2 <- xxx2[keep1 == TRUE & keep2 == TRUE, ]
	if(is.null(dim(xxx2))) { xxx2 <- as.matrix(xxx2) }
	
	# count the total number of included genotyped sites at this point
	n_sites <- nrow(xxx1)
	
	# combine the two matrices to find sites that are variant w/in and between the two pops
	xxx_combined <- cbind(xxx1, xxx2)
	variant_sites <- apply(xxx_combined, 1, polymorphic_function2)
	
	# keep only variant sites
	xxx1_variant <- xxx1[variant_sites > 1, ]
	if(is.null(dim(xxx1_variant))) { xxx1_variant <- as.matrix(xxx1_variant) }
	xxx2_variant <- xxx2[variant_sites > 1, ]
	if(is.null(dim(xxx2_variant))) { xxx2_variant <- as.matrix(xxx2_variant) }
	
	# count the number of variant sites
	n_variant_sites <- nrow(xxx1_variant)
	
	# loop for each polymorphic site to calculate dxy
	dxy_all <- list()
	for(a in 1:nrow(xxx1_variant)) {
		a_rep1 <- as.character(xxx1_variant[a,])
		a_rep2 <- as.character(xxx2_variant[a,])
		
		# remove missing
		a_rep1 <- a_rep1[a_rep1 != "./."]
		a_rep2 <- a_rep2[a_rep2 != "./."]
		
		# measure proportion of reference allele 
		a_ref1 <- (length(a_rep1[a_rep1 == "0/0"]) * 2 + length(a_rep1[a_rep1 == "0/1"]) * 1) / (length(a_rep1) * 2)
		a_ref2 <- (length(a_rep2[a_rep2 == "0/0"]) * 2 + length(a_rep2[a_rep2 == "0/1"]) * 1) / (length(a_rep2) * 2)
		
		# calc dxy
		dxy_all[[a]] <- a_ref1 * (1 - a_ref2) + a_ref2 * (1 - a_ref1)
	}
	dxy_all <- sum(unlist(dxy_all)) / n_sites
	
	
	# loop for each polymorphic site to calculate fst
	numerator_all <- list()
	denominator_all <- list()
	for(a in 1:nrow(xxx1_variant)) {
		a_rep1 <- as.character(xxx1_variant[a,])
		a_rep2 <- as.character(xxx2_variant[a,])
		
		# remove missing
		a_rep1 <- a_rep1[a_rep1 != "./."]
		a_rep2 <- a_rep2[a_rep2 != "./."]

		# number of individuals per population
		pop1_ind_count <- length(a_rep1) 
		pop2_ind_count <- length(a_rep2)
		
		# non-reference allele counts
		alt_allele_count1 <- (2 * length(a_rep1[a_rep1 == "1/1"]) + 1 * length(a_rep1[a_rep1 == "0/1"]))
		alt_allele_count2 <- (2 * length(a_rep2[a_rep2 == "1/1"]) + 1 * length(a_rep2[a_rep2 == "0/1"]))

		# total allele counts
		all_allele_count1 <- 2 * length(a_rep1)
		all_allele_count2 <- 2 * length(a_rep2)

		# expected heterozygosity for each population
		expected_het1 <- (alt_allele_count1 * (all_allele_count1 - alt_allele_count1)) / 
			(all_allele_count1 * (all_allele_count1 - 1))
		expected_het2 <- (alt_allele_count2 * (all_allele_count2 - alt_allele_count2)) / 
			(all_allele_count2 * (all_allele_count2 - 1))

		# find the fst numerator and denominator values for this snp (they all get summed and divided for 
		# the final estimate)
		numerator_all[[a]] <- (alt_allele_count1 / (2 * pop1_ind_count) - 
			alt_allele_count2 / (2 * pop2_ind_count))^2 - (expected_het1 / (2 * pop1_ind_count)) - 
			(expected_het2 / (2 * pop2_ind_count))
		denominator_all[[a]] <- numerator_all[[a]] + expected_het1 + expected_het2		
	}
	# calculate total fst for this window
	fst_all <- sum(unlist(numerator_all)) / sum(unlist(denominator_all))
	
	# write to output for dxy and fst
	output_rep1 <- c(popname1, popname2, "Dxy", strsplit(filename, ":")[[1]][1], 
				as.numeric(strsplit(strsplit(filename, ":")[[1]][2], "-")[[1]][1]),
				as.numeric(strsplit(strsplit(filename, "-")[[1]][2], ".simple")[[1]][1]),
				n_sites, n_variant_sites, dxy_all)
	output_rep2 <- c(popname1, popname2, "Fst", strsplit(filename, ":")[[1]][1], 
				as.numeric(strsplit(strsplit(filename, ":")[[1]][2], "-")[[1]][1]),
				as.numeric(strsplit(strsplit(filename, "-")[[1]][2], ".simple")[[1]][1]),
				n_sites, n_variant_sites, fst_all)
	write(output_rep1, file=outname, append=T, ncolumns=9, sep="\t")
	write(output_rep2, file=outname, append=T, ncolumns=9, sep="\t")
	
}


		