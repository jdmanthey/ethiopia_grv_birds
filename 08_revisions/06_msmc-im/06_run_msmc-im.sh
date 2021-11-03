
cd /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input

source activate py39

cd cossypha

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.23e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/cossypha/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt

cd ../crithagra

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.062e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/crithagra/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt

cd ../melaenornis

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.151e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/melaenornis/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt

cd ../parophasma

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.074e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/parophasma/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt

cd ../turdus

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.123e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/turdus/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt

cd ../zosterops

MSMC_IM.py -beta 1e-8,1e-6 -mu 2.294e-9 -o /lustre/scratch/jmanthey/03_ethiopia_popgen/07_msmc_input/zosterops/output.msmcim -p 1*2+20*1+1*2+1*3 --printfittingdetails --plotfittingdetails --xlog msmc.final.txt



