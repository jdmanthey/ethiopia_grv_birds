# maximum clade credibility tree (simple summarization) from all gene trees using dendropy:
# gives info about which nodes have support from what proportion of gene trees
sumtrees.py --output=ethiopia_100kbp.tre --min-clade-freq=0.05 ethiopia.trees 
