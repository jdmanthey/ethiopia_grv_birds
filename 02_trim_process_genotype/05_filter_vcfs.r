# use this r script to make a shell filtering script for each file in the vcf list
	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/02_ethiopia_popgen"
	directory_name <- "erv_scripts"
	queue <- "omni"
	cluster <- "quanah"
	vcftools_options2 <- "--max-missing 0.005 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --remove-indels"
	vcftools_options3 <- "--max-missing 1 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --remove-indels"
	
	# directory names and genus names of species (make sure directories are made first)
	genera <- c("cossypha", "melaenornis", "parophasma", "serinus", "turdus", "zosterops")
	
	# read in the vcf list file
	vcf_list <- "vcf_list.txt"
	vcf <- scan("vcf_list.txt", what="character")
	# define output base name
	out_name <- sapply(strsplit(vcf, "\\.g\\."), "[[", 1)

	# make directories
	dir.create(directory_name)
	
	# write the two helper files
	write(vcf, file=paste(directory_name, "/helper7.txt", sep=""), ncolumns=1)
	write(out_name, file=paste(directory_name, "/helper8.txt", sep=""), ncolumns=1)

	# write the array script
	a.script <- paste(directory_name, "/filter_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "filter_vcf", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 1", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=12:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=8G", file=a.script, append=T)
	write(paste("#$ -t 1:", length(vcf), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("input_array=$( head -n${SGE_TASK_ID} helper7.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("output_array=$( head -n${SGE_TASK_ID} helper8.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	vcftools_total <- paste("vcftools --vcf ", project_directory, "/03_vcf/${input_array} ", vcftools_options2, 
		" --recode --recode-INFO-all --out ", project_directory, "/03_vcf/${output_array}", sep="")
	write(vcftools_total, file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("bgzip ", project_directory, "/03_vcf/${output_array}.recode.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("tabix ", project_directory, "/03_vcf/${output_array}.recode.vcf.gz", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	# filter twice for each genus
	for(a in 1:length(genera)) {
		write("", file=a.script, append=T)
		vcftools_total <- paste("vcftools --vcf ", project_directory, "/03_vcf/${input_array} ", 
			vcftools_options2, 
			" --keep ", genera[a], ".txt --recode --recode-INFO-all --out ", project_directory,
			"/03_vcf/", genera[a], "/${output_array}__msmc", sep="")
		write(vcftools_total, file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("bcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n \' ", project_directory, 
			"/03_vcf/", genera[a], "/${output_array}__msmc.recode.vcf > ", project_directory, 
			"/03_vcf/", genera[a], "/${output_array}__msmc.simple.vcf", sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("rm ", project_directory, "/03_vcf/", genera[a], "/${output_array}__msmc.recode.vcf",
			sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		vcftools_total <- paste("vcftools --vcf ", project_directory, "/03_vcf/${input_array} ", 
			vcftools_options3, 
			" --keep ", genera[a], ".txt --recode --recode-INFO-all --out ", project_directory, 
			"/03_vcf/", genera[a], "/${output_array}__struc", sep="")
		write(vcftools_total, file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("bcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n \' ", project_directory, 
			"/03_vcf/", genera[a], "/${output_array}__struc.recode.vcf > ", project_directory, 
			"/03_vcf/", genera[a], "/${output_array}__struc.simple.vcf", sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("rm ", project_directory, "/03_vcf/", genera[a], "/${output_array}__struc.recode.vcf",
			sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
	}
			
