# make a gff file for each of the cds features in the overall gff file

x <- read.table("GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.gff", stringsAsFactors=F, sep="\t", fill=T)

# keep only cds
x2 <- x[x[,3] == "CDS",]

# change row names
rownames(x2) <- seq(from=1, to=nrow(x2), by=1)

# keep only basic attribute data
x2[,9] <- sapply(strsplit(x2[,9], ";"), "[[", 1)


# write the cds to file
write.table(x2, file="finch_cds2.gff", sep="\t", row.names=F, col.names=F, quote=F)
