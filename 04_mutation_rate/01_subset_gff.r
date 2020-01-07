# make a gff file for each of the cds features in the overall gff file

x <- read.table("GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.gff", stringsAsFactors=F, sep="", fill=T)

# keep only cds
x <- x[x[,3] == "CDS",]

# identify the unique genes
gene_list <- sapply(strsplit(sapply(strsplit(x[,9], ";"), "[[", 1), "="), "[[", 2)
unique_gene_list <- unique(gene_list)

# write the cds to file
write.table(x, file="finch_cds.gff", sep="\t", row.names=F, col.names=F, quote=F)
