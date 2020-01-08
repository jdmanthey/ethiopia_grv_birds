# remove some strange <NON_REF> calls as genotypes for some sites in the vcf file
grep -v "NON_REF" total.g.vcf > total2.g.vcf

# intersect the genes gff with the total vcf
bedtools intersect -u -a total2.g.vcf -b finch_cds2.gff > genes.vcf

# extract the header
grep "#" --line-buffered total2.g.vcf > header.vcf

# cat the bedtools extracted genes vcf file with the header
cat header.vcf genes.vcf > genes2.vcf

# filter the vcf with bedtools
vcftools --vcf genes2.vcf --max-missing 0.3 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 \
--remove-indels --recode --recode-INFO-all --out genes_filtered

# simplify the filtered vcf with bcftools
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' genes_filtered.recode.vcf > genes_filtered.simple.vcf
