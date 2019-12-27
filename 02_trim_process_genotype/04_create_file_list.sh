# move to directory with all the group genotype files and write a list of the total number of vcf files for use with 
# the next R script

ls -1 *vcf > vcf_list.txt

# make a vcf file that contains just the header of the vcf files (all are in common)

grep '#' NC_011462.1__a.g.vcf > header.vcf
