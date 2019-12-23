samtools depth -a 1_final.bam 2_final.bam 3_final.bam 4_final.bam 5_final.bam 6_final.bam 7_final.bam 8_final.bam \
9_final.bam 10_final.bam 11_final.bam 12_final.bam 13_final.bam 14_final.bam 15_final.bam 16_final.bam 17_final.bam \
18_final.bam 19_final.bam 20_final.bam 21_final.bam 22_final.bam 23_final.bam 24_final.bam 25_final.bam \
26_final.bam 27_final.bam 28_final.bam 29_final.bam 30_final.bam > \
eth_grv_coverage.txt

# break up the depth files into single column files for each individual (locations dropped)

while read -r name1 number1; do
	number2=$((number1 + 2));
  cut eth_grv_coverage.txt -f $number2 > ${name1}_depth.txt;
done < ethiopia_grv_popmap.txt

