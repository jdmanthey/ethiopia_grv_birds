# makes an array to convert all the vcf files into windowed simple vcfs that can be used to calculate stats
# use the helper5 and helper6 text files from the genotyping third step to help subset the dataset into windows
# make sure a "windows" directory is made in the 03_vcf directory in the project directory

	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/02_ethiopia_popgen"
	directory_name <- "ethiopia_windows"
	queue <- "omni"
	cluster <- "quanah"
	
	# define window size
	window_size <- 500000
	
	# interval per file definitions
	helper5 <- scan("helper5.txt", what="character")
	# file name definitions
	helper6 <- scan("helper6.txt", what="character")
	helper6 <- paste(helper6, ".recode.vcf.gz", sep="")
	
	# make directories
	dir.create(directory_name)
	
	# define intervals and write to helper files
	helpers <- c()
	for(a in 1:length(helper6)) {
		a_start <- as.numeric(strsplit(strsplit(helper5[a], ":")[[1]][2], "-")[[1]][1])
		a_end <- a_start + window_size - 1
		a_max <- as.numeric(strsplit(helper5[a], "-")[[1]][2])
		a_windows <- ceiling((a_max - a_start) / window_size)
		a_chromosome <- strsplit(helper5[a], ":")[[1]][1]
		
		# loop for defining helper info for each window
		for(b in 1:a_windows) {
			if(b == a_windows) {
				a_end <- a_max
			}
			b_helper <- c(helper6[a], a_chromosome, a_start, a_end, paste(a_chromosome, ":", a_start, "-", a_end, sep=""))
			helpers <- rbind(helpers, b_helper)
			a_start <- a_start + window_size
			a_end <- a_end + window_size
		}
	}
	write(helpers[,1], file=paste(directory_name, "/helper9.txt", sep=""), ncolumns=1)
	write(helpers[,5], file=paste(directory_name, "/helper10.txt", sep=""), ncolumns=1)
	
	
	
	# write the array script
	a.script <- paste(directory_name, "/window_split_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "window_split", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 1", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=01:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=8G", file=a.script, append=T)
	write(paste("#$ -t 1:", nrow(helpers), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("input_array=$( head -n${SGE_TASK_ID} helper9.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("index_array=$( head -n${SGE_TASK_ID} helper10.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	tabix_command <- paste("tabix ", project_directory, "/03_vcf/${input_array} ${index_array} > ", 
		"temp${SGE_TASK_ID}.vcf", sep="")
	write(tabix_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	cat_command <- paste("cat ", project_directory, "/03_vcf/header.vcf ", 
		"temp${SGE_TASK_ID}.vcf > ", project_directory, "/03_vcf/windows/${index_array}.vcf", sep="")
	write(cat_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("rm temp${SGE_TASK_ID}.vcf"), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("bcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n\' ", 
		project_directory, "/03_vcf/windows/${index_array}.vcf > ", 
		project_directory, "/03_vcf/windows/${index_array}.simple.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("rm ", project_directory, "/03_vcf/windows/${index_array}.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)

		
	
	
	
