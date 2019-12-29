# run these commands in the windows directory with the output
# header for each combined file
grep 'pop1' NC_011462.1:1-500000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' NC_011462.1:1-500000__stats.txt > ../window_fst.txt
grep 'pop1' NC_011462.1:1-500000__stats.txt > ../window_dxy.txt

# add the relevant statistics to each file
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *txt ); do grep 'Dxy' $i >> ../window_Dxy.txt; done
for i in $( ls *txt ); do grep 'Fst' $i >> ../window_Fst.txt; done
