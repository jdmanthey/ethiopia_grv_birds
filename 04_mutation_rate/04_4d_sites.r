
library(stats)
library(Biostrings)
library(seqinr)
library(rphast)

# make the directory for filtering the fasta files
dir.create("_filtered")

# make the 4 fold degenerate sites output directory
dir.create("_4d_output")

# list all input files
x_files <- list.files("fasta", full.names=T)

# loop to filter out any codons with any Ns 
for(a in 1:length(x_files)) {
	a_rep <- read.msa(x_files[a])
	# make sure the alignment is divisible by 3
	if(nchar(a_rep[[1]][1]) %% 3 == 0) {
		# loop for each individual
		missing <- c()
		for(b in 1:length(a_rep[[1]])) {
			missing <- c(missing, grep("N", strsplit(a_rep[[1]][b], "")[[1]]))
		}
		# remove redundant sites
		missing <- sort(unique(missing))
		
		# keep going if missing has any sites
		if(!is.null(missing)) {
			missing_codons <- unique(ceiling(missing / 3))
			# define the missing sites for these codons
			missing_sites <- sort(c(missing_codons * 3 - 2, missing_codons * 3 - 1, missing_codons * 3))
			# define the sites to keep
			sites_to_keep <- seq(from=1, to=nchar(a_rep[[1]][1]), by=1)[seq(from=1, to=nchar(a_rep[[1]][1]), by=1) %in% missing_sites == FALSE]
			
		} else {
			# keep all sites if none missing
			sites_to_keep <- seq(from=1, to=nchar(a_rep[[1]][1]), by=1)
		}
		
		# write output
		a_name <- paste("_filtered/", sapply(strsplit(x_files[a], "/"), "[[", 2), sep="")
		# loop for each individual
		for(b in 1:length(a_rep[[1]])) {
			if(b == 1) {
				write(paste(">", a_rep[[2]][b], sep=""), file=a_name, ncolumns=1)
			} else {
				write(paste(">", a_rep[[2]][b], sep=""), file=a_name, ncolumns=1, append=T)
			}
			# filter output sites and write
			output_sites <- paste(strsplit(a_rep[[1]][b], "")[[1]][sites_to_keep], collapse="")
			write(output_sites, file=a_name, ncolumns=1, append=T)
		}
	}
}



# list all filtered files
x_files <- list.files("_filtered", full.names=T)

# loop to
# read in multiple sequence alignments to get the four-fold degenerate sites
# for a later phylogeny
for(a in 1:length(x_files)) {
	a_rep <- read.msa(x_files[a])
	
	# make sure the alignment is divisible by 3
	if(nchar(a_rep[[1]][1]) %% 3 == 0) {
		# define the coding feature in the alignment (they are all trimmed to be in 1st position starting the codons)
		a_feature <- feat(seqname="Taeniopygia_guttata", feature="CDS", start=1, end=ncol(a_rep))
		# get the 4 fold degenerate sites
		a_4d_rep <- get4d.msa(a_rep, a_feature)
		
		# determine basename of output
		a_name <- strsplit(sapply(strsplit(x_files[a], "/"), "[[", 2), ".fasta")[[1]]
		# write the alignment
		write.msa(a_4d_rep, file=paste("_4d_output/", a_name, "_4d.fasta", sep=""), format="FASTA")
	}
}



# list all the 4d alignments output
x_files <- list.files("_4d_output", full.names=T)
# loop to read in all alignments and concatenate
cossypha <- list()
melaenornis <- list()
parophasma <- list()
serinus <- list()
turdus <- list()
zosterops <- list()
taeniopygia <- list()

for(a in 1:length(x_files)) {
	a_rep <- readDNAStringSet(x_files[a])
	cossypha[[a]] <- as.character(a_rep)[1]
	melaenornis[[a]] <- as.character(a_rep)[2]
	parophasma[[a]] <- as.character(a_rep)[3]
	serinus[[a]] <- as.character(a_rep)[4]
	turdus[[a]] <- as.character(a_rep)[5]
	zosterops[[a]] <- as.character(a_rep)[6]
	taeniopygia[[a]] <- as.character(a_rep)[7]
}
cossypha <- paste(unlist(cossypha), collapse="")
melaenornis <- paste(unlist(melaenornis), collapse="")
parophasma <- paste(unlist(parophasma), collapse="")
serinus <- paste(unlist(serinus), collapse="")
turdus <- paste(unlist(turdus), collapse="")
zosterops <- paste(unlist(zosterops), collapse="")
taeniopygia <- paste(unlist(taeniopygia), collapse="")

output_name <- "_total_4d_sites.fasta"
write(">cossypha", file=output_name, ncolumns=1)
write(cossypha, file=output_name, ncolumns=1, append=T)
write(">melaenornis", file=output_name, ncolumns=1, append=T)
write(melaenornis, file=output_name, ncolumns=1, append=T)
write(">parophasma", file=output_name, ncolumns=1, append=T)
write(parophasma, file=output_name, ncolumns=1, append=T)
write(">serinus", file=output_name, ncolumns=1, append=T)
write(serinus, file=output_name, ncolumns=1, append=T)
write(">turdus", file=output_name, ncolumns=1, append=T)
write(turdus, file=output_name, ncolumns=1, append=T)
write(">zosterops", file=output_name, ncolumns=1, append=T)
write(zosterops, file=output_name, ncolumns=1, append=T)
write(">taeniopygia", file=output_name, ncolumns=1, append=T)
write(taeniopygia, file=output_name, ncolumns=1, append=T)
