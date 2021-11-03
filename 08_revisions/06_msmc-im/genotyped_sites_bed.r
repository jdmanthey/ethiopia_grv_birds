# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# read in input file
input_file <- read.table(args[1], sep="\t", stringsAsFactors=F)

# define output file
output_name <- paste0(substr(args[1], 1, nchar(args[1]) - 16), ".coverage.mask")

# no scientific notation
options(scipen=999)

# determine unique chromosomes
chr <- unique(input_file[,1])

# loop through each individual chromosome
for(a in 1:length(chr)) {
        # subset this chromosome
        a_rep <- input_file[input_file[,1] == chr[a],2]
        
        # get the differences between each of the numbers (and add placeholder 1 at beginning)
        a_diff <- c(1, diff(a_rep))
        
        # create numbered index of a_diff object
        a_index <- seq(from=1, to=length(a_diff), by=1)[a_diff != 1]
        
        # get the first coordinate for each interval
        a_coord1 <- a_rep[c(1, a_index)]
        
        # get the second coordinate for each interval
        a_coord2 <- a_rep[c(a_index - 1, length(a_rep))]
        
        # subtract one from the first coordinate for each interval to make zero-based
        a_coord1 <- a_coord1 - 1
        
        # make output table
        output <- data.frame(chrom=as.character(chr[a]), coord1=as.numeric(a_coord1), coord2=as.numeric(a_coord2))
        
        # write output table
        if(a == 1) {
                write.table(output, file=output_name, sep="\t", quote=F, col.names=F, row.names=F)
        } else {
                write.table(output, file=output_name, sep="\t", quote=F, col.names=F, row.names=F, append=T)
        }
        
}       


