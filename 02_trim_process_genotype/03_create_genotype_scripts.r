# popmap = individual base names of fastq files, one line per individual

# make sure reference is indexed with bwa and samtools before use, and use CreateSequenceDictionary in GATK 
	options(scipen=999)

	project_directory <- "/lustre/scratch/jmanthey/02_ethiopia_popgen"
	directory_name <- "eth_grv_scripts"
	reference_genome_location <- "/home/jmanthey/references/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna"
	queue <- "omni"
	cluster <- "quanah"
	output_name <- "eth_grv"
	popmap <- "ethiopia_grv_popmap.txt"
	individuals <- read.table(popmap, sep="\t")
	faidx <- read.table("GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna.fai", sep="\t", stringsAsFactors=F)
	
	min_scaffold_size <- 2000000
	max_genotype_job_size <- 10000000

	# make directories
	dir.create(directory_name)
	dir.create(paste(directory_name, "/01_gatk_split", sep=""))
	dir.create(paste(directory_name, "/02b_gatk_database", sep=""))
	dir.create(paste(directory_name, "/03b_group_genotype_database", sep=""))

	# subset the index
	faidx <- faidx[faidx[,2] >= min_scaffold_size, ]

	
	
	
	# write the helper files for the 1st genotyping step
	for(a in 1:nrow(individuals)) {
		if(a == 1) {
			helper1 <- cbind(faidx[,1], rep(a, nrow(faidx)))
		} else {
			helper1 <- rbind(helper1, cbind(faidx[,1], rep(a, nrow(faidx))))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper1[,1], file=paste(directory_name, "/01_gatk_split", "/helper1.txt", sep=""), ncolumns=1)
	# write the individual to genotype for each job
	write(helper1[,2], file=paste(directory_name, "/01_gatk_split", "/helper2.txt", sep=""), ncolumns=1)
	
	
	
	# step 1
	# genotype all individuals using GATK, one array job using the above two helper files

	a.script <- paste(directory_name, "/01_gatk_split/step1_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "step1", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 8", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=48:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=8G", file=a.script, append=T)
	write(paste("#$ -t 1:", nrow(helper1), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SGE_TASK_ID} helper1.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("ind_array=$( head -n${SGE_TASK_ID} helper2.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	#gatk 4.0
	a_name <- paste(project_directory, "/01_bam_files/", "${ind_array}", "_final.bam", sep="")
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx64g" HaplotypeCaller -R ', reference_genome_location, " -I ", a_name, " -ERC GVCF -O ", project_directory, "/02_vcf/", "${chr_array}", "._${ind_array}_.g.vcf", " --QUIET --intervals ", "${chr_array}", sep="")
	write(gatk_command, file=a.script, append=T)
	
	
	
	
	# write the helper files for the 2nd genotyping step
	# write the chromosome/scaffold to database for each job
	write(faidx[,1], file=paste(directory_name, "/02b_gatk_database", "/helper3.txt", sep=""), ncolumns=1)

	
	
	# step 2
	# create genotyping database for each of the chromosomes
	a.script <- paste(directory_name, "/02b_gatk_database/step2_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "step2", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 6", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=48:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=10G", file=a.script, append=T)
	write(paste("#$ -t 1:", nrow(faidx), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SGE_TASK_ID} helper3.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
		
	#make list of all vcfs to database
	for(b in 1:nrow(individuals)) {
		if(b == 1) {
			vcf_total <- paste("-V ", project_directory, "/02_vcf/", "${chr_array}", "._", b, "_.g.vcf", sep="")
		} else {
			vcf_total <- paste(vcf_total, " -V ", project_directory, "/02_vcf/", "${chr_array}", "._", b, "_.g.vcf", sep="")
		}
	}
	
	#gatk 4.0
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx60g" GenomicsDBImport ', vcf_total, " --genomicsdb-workspace-path ", project_directory, "/02_vcf/", "${chr_array}", " -L ", "${chr_array}", sep="")
				
	write(gatk_command, file=a.script, append=T)
		




	# write the helper files for the 3rd genotyping step
	helper <- c()
	job_suffixes <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
	for(a in 1:nrow(faidx)) {
		a_rep <- faidx[a,]
		job_number <- ceiling(a_rep[1,2] / max_genotype_job_size)
		job_number <- job_suffixes[1:job_number]
		
		if(length(job_number) > 1) {
			# define interval to group genotype if more than one job
			interval_start <- 1
			interval_end <- max_genotype_job_size
			for(b in 1:length(job_number)) {
				helper <- rbind(helper, c(faidx[a,1], 
					paste(faidx[a,1], ":", interval_start, "-", interval_end, sep=""),
					paste(faidx[a,1], "__", job_number[b], sep="")
					))
					interval_start <- interval_start + max_genotype_job_size
				if((b+1) != length(job_number)) {
					interval_end <- interval_end + max_genotype_job_size
				} else {
					interval_end <- faidx[a,2]
				}
			}
		} else {
			helper <- rbind(helper, c(faidx[a,1], 
				paste(faidx[a,1], ":", 1, "-", faidx[a,2], sep=""),
				paste(faidx[a,1], sep="")
				))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper[,1], file=paste(directory_name, "/03b_group_genotype_database", "/helper4.txt", sep=""), ncolumns=1)
	# write the interval range to genotype for each job
	write(helper[,2], file=paste(directory_name, "/03b_group_genotype_database", "/helper5.txt", sep=""), ncolumns=1)
	# write the output base name for each job
	write(helper[,3], file=paste(directory_name, "/03b_group_genotype_database", "/helper6.txt", sep=""), ncolumns=1)
	
	
	


	

	# step 3
	# group genotyping for each interval
	a.script <- paste(directory_name, "/03b_group_genotype_database/step3_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "step3", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 10", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=48:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=8G", file=a.script, append=T)
	write(paste("#$ -t 1:", nrow(helper), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SGE_TASK_ID} helper4.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("interval_array=$( head -n${SGE_TASK_ID} helper5.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SGE_TASK_ID} helper6.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)

	#gatk 4.0
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx80g" GenotypeGVCFs -R ', reference_genome_location, " -V gendb://", project_directory, "/02_vcf/", "${chr_array}", " --include-non-variant-sites -O ", project_directory, "/03_vcf/", "${name_array}", ".g.vcf", " -L ", "${interval_array}", sep="")
	write(gatk_command, file=a.script, append=T)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
