# list all input files in popmap
ind_list <- scan("popmap.txt", what="character")
species <- unique(paste(sapply(strsplit(ind_list, "_"), "[[", 1), sapply(strsplit(ind_list, "_"), "[[", 2), sep="_"))
species_all <- paste(sapply(strsplit(ind_list, "_"), "[[", 1), sapply(strsplit(ind_list, "_"), "[[", 2), sep="_")


# define parameters
# species order = cossypha, crithagra, melaenornis, parophasma, turdus, zosterops
# gen time = 2 * age of sexual maturity in closely-related species from genomics.senescence.info/species
# age of sex. maturity for each = one year
# method citation = doi: 10.1016/j.cub.2015.03.047
gen <- 2
mu <- c(2.23024e-9, 2.061943e-9, 2.151193e-9, 2.074098e-9, 2.122962e-9, 2.294e-9)* gen # needs to be per generation
# species number (number per each species)
# 5, 5, 6, 3, 6, 5
mu <- c(rep(mu[1], 5), rep(mu[2], 5), rep(mu[3], 6), rep(mu[4], 3), rep(mu[5], 6), rep(mu[6], 5))
min_age <- 000
max_age <- 500000
plotting_age_range <- 1000 # years for x axis


# loop through each output directory
output <- list()
for(a in 1:length(ind_list)) {

	# identify main output file
	a_main <- paste("output/", ind_list[a], ".final.txt", sep="")

	
	# read in main file
	a_rep <- read.table(a_main, sep="\t", header=T)
	# rearrange main file for plotting lines
	for(d in 1:nrow(a_rep)) {
		if(d == 1) {
			a_rep2 <- rbind(c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		} else {
			a_rep2 <- rbind(a_rep2, c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		}
	}
	a_rep <- a_rep2
	# scale by mutation rate
	a_rep[,1] <- a_rep[,1] / mu[a] * gen
	a_rep[,2] <- (1 / a_rep[,2]) / (2 * mu[a])
	# remove very young and old time frames prone to error
	a_rep <- a_rep[a_rep[,1] >= min_age & a_rep[,1] <= max_age,]
	# scale by plotting age range and pop size range
	a_rep <- a_rep / plotting_age_range
	# add to output list
	output[[a]] <- a_rep
	
	
}

# define current pop sizes 
current_pop <- c()
for(a in 1:length(output)) {
	current_pop <- c(current_pop, output[[a]][1,2])
}

# harmonic mean of pop sizes from most recent to 200k years ago
harmonic_pop <- c()
for(a in 1:length(output)) {
	out_rep <- output[[a]]
	
	# define time series
	time_series <- seq(from=as.integer(out_rep[2,1])+1, to=200, by=1)
	# time series pops
	time_pops <- c()
	for(b in 1:length(time_series)) {
		time_pops <- c(time_pops, out_rep[time_series[b] < out_rep[,1],][1,2])
	}
	# harmonic mean of this individual
	harm_rep <- length(time_pops) / sum((1 / time_pops))
	
	# add to output element
	harmonic_pop <- c(harmonic_pop, harm_rep)
}

# define output
output <- data.frame(individual=as.character(ind_list), current_pop=as.numeric(current_pop), harmonic_pop=as.numeric(harmonic_pop))

#write output
write.table(output, file="pop_sizes.txt", sep="\t", row.names=F, quote=F)
