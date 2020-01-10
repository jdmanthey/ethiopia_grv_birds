# for a comparison tree I pruned the phylogeny of all passeriformes bird families (Earth history and the passerine superradiation)
# here i pruned to six taxa that were the same as the genes used in this study or representatives from their 
# respective families to include
# here, cossypha and melaenornis used the same mutation rates


# Cossypha = Muscicapa
# Crithagra = Fringilla
# Melaenornis = Muscicapa
# Parophasma = Curruca
# Turdus = Turdus
# Zosterops = Zosterops
# Taeniopygia = Estrilda

# used the following R code to extract taxa of interest (or taxa from same family)

library(ape)

x <- read.nexus(file="Figure-1.txt")

tips_to_keep <- c("Muscicapa_striata", "Fringilla_montifringilla", "Curruca_nana", "Turdus_albicollis", "Zosterops_everetti", "Estrilda_melpoda")

x2 <- keep.tip(x, tips_to_keep)

write.tree(x2, file="pruned_tree.tre")
