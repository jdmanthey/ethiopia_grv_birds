cd /lustre/scratch/jmanthey/02_ethiopia_popgen/03_vcf

# make a new directory
mkdir extract_genes

# put the vcf header in a new file
grep '#' NC_011462.1__a.g.vcf >> extract_genes/total.g.vcf

# add the vcf contents all to one file for bedtools to use
for i in $( ls *g.vcf ); do
	grep -v '#' $i >> extract_genes/total.g.vcf
done

