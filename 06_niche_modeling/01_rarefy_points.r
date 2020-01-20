# rarefy unique function rarefies your file to only unique occurrence points
# input for rarefy_unique is a csv with columns lat, long
# example usage:
# source("rarefy_update.r")
# rarefy_unique("test.csv")

rarefy_unique <- function(file) {
	x <- as.matrix(read.csv(file))
	
	# keep only unique coordinates
	x <- rbind(c("Lat","Long"), unique(x))
	
	# write output
	write.table(x, file = paste(substr(file, 1, nchar(file) - 4), "_unique.csv", sep=""), row.names=F,
		col.names=F, quote=F, sep=",")
}


# rarefy function rarefies your file to points at least distance.km kilometers apart
# input for rarefy is a csv with columns lat, long
# example usage:
# source("rarefy_update.r")
# rarefy("test.csv", 5)

rarefy <- function(file, distance.km) {
	require(fossil)
	x <- as.matrix(read.csv(file))
	
	# ensure distance is numeric
	distance.km <- as.numeric(distance.km)
	
	# keep only unique coordinates
	x <- unique(x)
	
	# counter for continuing rarefication
	reps <- nrow(x) - 1
	
	# create output file
	output <- c()
	
	# create initial object for testing
	x2 <- x
	
	
	# if large dataset, subset first to narrow down
	if(reps > 1000) {
		writeLines(paste("There are many points, subsetting"))

		x3 <- x2
		s.reps <- as.integer(reps / 1000) + 1
		
		# create 1,000 point subsets
		subsets <- list()
		for(a in 1:s.reps) {
			if(a == s.reps) {
				subsets[[a]] <- x3
			} else {
				subsets[[a]] <- x3[1:1000, ]
				x3 <- x3[1001:nrow(x3),]
			}
			
		}
		
		# distance filter each subset
		output.temp <- list()
		for(a in 1:s.reps) {
			writeLines(paste("Subset a", a*1000, "of", s.reps*1000))
			reps.temp <- nrow(subsets[[a]])
			out.temp <- c()
			a.rep.temp <- subsets[[a]]
			# distance functions while > 20 points left (arbitrary number)
			while(reps.temp > 20) {
				a.rep <- as.vector(a.rep.temp[1,])
				a.rep.temp <- a.rep.temp[2:nrow(a.rep.temp),]
				out.temp <- rbind(out.temp, a.rep)
				testing <- apply(a.rep.temp, 1, distance, yy.value = a.rep)
				a.rep.temp <- a.rep.temp[testing > distance.km, ]
				if(!is.null(nrow(a.rep.temp))) {
					reps.temp <- nrow(a.rep.temp) - 1
				} else {
					reps.temp <- 1
				}
				
			}
			out.temp <- rbind(out.temp, a.rep.temp)
			output.temp[[a]] <- out.temp
		}
		# sum total number of points left 
		total <- 0
		for(a in 1:s.reps) {
			total <- total + nrow(output.temp[[a]])
		}
		writeLines(paste("There are still", total, "points left"))
		
		# take all output and put in a new matrix
		new.x2 <- c()
		for(a in 1:s.reps) {
			new.x2 <- rbind(new.x2, output.temp[[a]])
		}
		x2 <- new.x2
		reps <- nrow(x2) - 1
	}
	
	# if still large dataset, subset again to narrow down
	if(reps > 1000) {
		writeLines(paste("There are many points, subsetting again"))

		x3 <- x2
		s.reps <- as.integer(reps / 1000) + 1
		
		# create 1,000 point subsets
		subsets <- list()
		for(a in 1:s.reps) {
			if(a == s.reps) {
				subsets[[a]] <- x3
			} else {
				subsets[[a]] <- x3[1:1000, ]
				x3 <- x3[1001:nrow(x3),]
			}
			
		}
		
		# distance filter each subset
		output.temp <- list()
		for(a in 1:s.reps) {
			writeLines(paste("Subset b", a*1000, "of", total))
			reps.temp <- nrow(subsets[[a]])
			out.temp <- c()
			a.rep.temp <- subsets[[a]]
			# distance functions while > 20 points left (arbitrary number)
			while(reps.temp > 20) {
				a.rep <- as.vector(a.rep.temp[1,])
				a.rep.temp <- a.rep.temp[2:nrow(a.rep.temp),]
				out.temp <- rbind(out.temp, a.rep)
				testing <- apply(a.rep.temp, 1, distance, yy.value = a.rep)
				a.rep.temp <- a.rep.temp[testing > distance.km, ]
				if(!is.null(nrow(a.rep.temp))) {
					reps.temp <- nrow(a.rep.temp) - 1
				} else {
					reps.temp <- 1
				}
			}
			out.temp <- rbind(out.temp, a.rep.temp)
			output.temp[[a]] <- out.temp
		}
		# sum total number of points left 
		total <- 0
		for(a in 1:s.reps) {
			total <- total + nrow(output.temp[[a]])
		}
		writeLines(paste("There are still", total, "points left"))
		
		# take all output and put in a new matrix
		new.x2 <- c()
		for(a in 1:s.reps) {
			new.x2 <- rbind(new.x2, output.temp[[a]])
		}
		x2 <- new.x2
		reps <- nrow(x2) - 1
	}
	
	
	# while loop to keep going until filtered through entire file
	while(reps >= 1) {
		if(reps %% 50 == 0) {
			writeLines(paste("Final Countdown:", reps))
		}
		# distance functions while > 2 points left
		if(reps > 1) {
			a.rep <- as.vector(x2[1,])
			x2 <- x2[2:nrow(x2),]
			output <- rbind(output, a.rep)
			testing <- apply(x2, 1, distance, yy.value = a.rep)
			x2 <- x2[testing > distance.km, ]
			if(!is.null(nrow(x2))) {
				reps <- nrow(x2) - 1
			} else {
				if(length(x2) > 0) {
					x2 <- rbind(as.vector(x2), as.vector(x2))
					reps <- 1
				} else {
					reps <- 0
				}
			}
		} else if(reps == 1) { # finish if there is only one point left
			x2 <- x2[2:nrow(x2),]
			output <- rbind(output, as.vector(x2))
			reps <- 0
		}
		
	}
	
	# write output
	output <- rbind(c("Lat", "Long"), output)
	write.table(output, file = paste(substr(file, 1, nchar(file) - 4), "_", distance.km, "km_rarefy.csv", 
		sep=""), row.names=F, col.names=F, quote=F, sep=",")
	
}


# function to be used with apply to get all distances from a matrix
# only used within other functions
distance <- function(xx.matrix, yy.value) {
	d.rep <- deg.dist(xx.matrix[2], xx.matrix[1], yy.value[2], yy.value[1])
	return(d.rep)
}
