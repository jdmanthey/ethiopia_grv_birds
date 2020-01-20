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
output_bootstraps <- list()
for(a in 1:length(ind_list)) {
	# identify bootstraps 
	a_bootstraps <- paste("output/", ind_list[a], "_b", seq(from=1, to=10, by=1), ".final.txt", sep="")
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
	
	# output for each bootstrap
	output_bootstraps[[a]] <- list()
	for(b in 1:length(a_bootstraps)) {
		b_rep <- read.table(a_bootstraps[b], sep="\t", header=T)
		# rearrange main file for plotting lines
		for(d in 1:nrow(b_rep)) {
			if(d == 1) {
				b_rep2 <- rbind(c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			} else {
				b_rep2 <- rbind(b_rep2, c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			}
		}
		b_rep <- b_rep2
		# scale by mutation rate
		b_rep[,1] <- b_rep[,1] / mu[a] * gen
		b_rep[,2] <- (1 / b_rep[,2]) / (2 * mu[a])
		# remove very young and old time frames prone to error
		b_rep <- b_rep[b_rep[,1] >= min_age & b_rep[,1] <= max_age,]
		# scale by plotting age range and pop size range
		b_rep <- b_rep / plotting_age_range
		# add to output list
		output_bootstraps[[a]][[b]] <- b_rep		
	}
}



# plot all separate
a_col <- "darkgreen"
par(mfrow=c(3,10))
par(mar=c(1,1,2,0.2))
for(a in 1:length(output)) {
	plot_name1 <- ind_list[a]
	plot(c(-1,1), xlim=c(20, 300), ylim=c(0,1000), pch=19, cex=0.01, log="x", xlab="", ylab="", main="", xaxt="n", yaxt="n")
	title(main=bquote(italic(.(plot_name1))), adj=0,line=0.5, cex.main=0.8)
	axis(side=2, at=c(0, 200, 400, 600, 800, 1000), labels=FALSE)		
	axis(side=1, at=c(10, 20, 50, 100, 200), labels=F)		
	
	
	# plot bootstraps
	for(b in 1:length(output_bootstraps[[1]])) {
		lines(output_bootstraps[[a]][[b]][,1], output_bootstraps[[a]][[b]][,2], col=a_col, lwd=0.3)
	}
	lines(output[[a]][,1], output[[a]][,2], col=a_col, lwd=3)
}




# plot just the main estimates for east and west of the GRV for each species in panels for figure 2
library(scales)
a_col <- c(alpha("darkgreen", 0.7), alpha("darkorange4", 0.7), alpha("darkgoldenrod", 0.7))
par(mfrow=c(6,2))
par(mar=c(1,1,2,0.2))
# reorder output2 for plotting all populations together in one panel
output2_pop <- c("cossypha_east", "cossypha_east", "cossypha_east", "cossypha_west", "cossypha_west",
					"crithagra_east", "crithagra_east", "crithagra_west", "crithagra_west", "crithagra_west",
					"mel_east", "mel_east", "mel_west","mel_west","mel_west","mel_east",
					"par_east", "par_east", "par_west",
					"turd_east", "turd_east", "turd_east", "turd_west", "turd_west", "turd_west",
					"zost_east", "zost_west", "zost_west", "zost_west", "zost_east")
pops <- c("cossypha_west", "cossypha_east", "crithagra_west", "crithagra_east", "mel_west", "mel_east", "par_west", "par_east", "turd_west", "turd_east", "zost_west", "zost_east")

# plot
for(a in 1:length(pops)) {
	pop_rep <- seq(from=1, to=length(output2_pop), by=1)[output2_pop == pops[a]]
	plot_name1 <- pops[a]
	plot(c(-1,1), xlim=c(20, 300), ylim=c(0,1000), pch=19, cex=0.01, log="x", xlab="", ylab="", main="", xaxt="n", yaxt="n")
	title(main=bquote(italic(.(plot_name1))), adj=0,line=0.5, cex.main=0.8)
	axis(side=2, at=c(0, 200, 400, 600, 800, 1000), labels=FALSE)		
	axis(side=1, at=c(20, 50, 100, 200), labels=F)	
	# loop for each plot in this panel
	for(b in 1:length(pop_rep)) {
		lines(output[[pop_rep[b]]][,1], output[[pop_rep[b]]][,2], col=a_col[b], lwd=2)
	}
	
}














