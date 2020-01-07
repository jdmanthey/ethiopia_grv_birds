# run the msmc utility script for each directory to create bootstrap replicates

for i in $( ls ); do 
cd $i;
~/multihetsep_bootstrap.py -n 10 -s 1000000 --chunks_per_chromosome 30 --nr_chromosomes 31 \
--seed 324324 bootstrap *txt;
cd ..;
done
