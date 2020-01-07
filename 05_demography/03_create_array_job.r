individuals <- c("C_semirufa_EB008","C_semirufa_EB009","C_semirufa_EB024","C_semirufa_EB044","C_semirufa_EB062","C_tristriata_EB007","C_tristriata_EB015","C_tristriata_EB043","C_tristriata_EB058","C_tristriata_EB064","M_chocolatinus_EB002","M_chocolatinus_EB003","M_chocolatinus_EB053","M_chocolatinus_EB059","M_chocolatinus_EB063","M_chocolatinus_EB067","P_galinieri_EB004","P_galinieri_EB017","P_galinieri_EB056","T_abyssinicus_EB001","T_abyssinicus_EB006","T_abyssinicus_EB010","T_abyssinicus_EB042","T_abyssinicus_EB045","T_abyssinicus_EB046","Z_poliogastrus_EB013","Z_poliogastrus_EB049","Z_poliogastrus_EB050","Z_poliogastrus_EB051","Z_poliogastrus_EB065")

# species
species <- c("C_semirufa", "C_tristriata", "M_chocolatinus", "P_galinieri", "T_abyssinicus", "Z_poliogastrus")

# half the heterozygosity (mean per genus) for the MSMC r parameter
heterozygosity <- c(0.00229627, 0.00132413, 0.00139479, 0.00043941, 0.00215914, 0.00147503)

directory <- "/lustre/scratch/jmanthey/02_ethiopia_popgen/04_demography"

options(scipen=999)

queue <- "omni"
cluster <- "quanah"

# determine all the input directories
directories <- c()
for(a in 1:length(individuals)) {
	directories <- c(directories, paste(directory, "/", individuals[a], "/*txt", sep=""))
	# bootstraps
	for(b in 1:10) {
		directories <- c(directories, paste(directory, "/", individuals[a], "/bootstrap_", b, "/*txt", sep=""))
	}
}

# determine the output names for each 
outputs <- substr(sapply(strsplit(directories, directory), "[[", 2), 2, nchar(sapply(strsplit(directories, directory), "[[", 2)) - 5)
outputs <- gsub("/bootstrap_", "_b", outputs)

# determine the heterozygosities to use for each 
heterozygosities <- heterozygosity[match(sapply(strsplit(outputs, "_EB"), "[[", 1), species)]

# write helper files
write(directories, file="helper1.txt", ncolumns=1)
write(outputs, file="helper2.txt", ncolumns=1)
write(heterozygosities, file="helper3.txt", ncolumns=1)

# write the array job
a_script <- "eth_grv_array.sh"
write("#!/bin/sh", file=a_script)
write("#$ -V", file=a_script, append=T)
write("#$ -cwd", file=a_script, append=T)
write("#$ -S /bin/bash", file=a_script, append=T)
write(paste("#$ -N ", "grv_dem", sep=""), file=a_script, append=T)
write(paste("#$ -q ", queue, sep=""), file=a_script, append=T)
write("#$ -pe sm 2", file=a_script, append=T)
write(paste("#$ -P ", cluster, sep=""), file=a_script, append=T)
write("#$ -l h_rt=48:00:00", file=a_script, append=T)
write("#$ -l h_vmem=8G", file=a_script, append=T)
write(paste("#$ -t 1:", length(directories), sep=""), file=a_script, append=T)
write("", file=a_script, append=T)

write("input_array=$( head -n${SGE_TASK_ID} helper1.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("output_array=$( head -n${SGE_TASK_ID} helper2.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("heter_array=$( head -n${SGE_TASK_ID} helper3.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)

msmc_command <- "~/msmc2_linux64bit -o ${output_array} -i 20 -t 2 -m ${heter_array} -p 1*2+20*1+1*2+1*3 ${input_array}"
write(msmc_command, file=a_script, append=T)



