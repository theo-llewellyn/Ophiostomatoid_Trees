library(tidyverse)
library(ape)

#drop tips in species for astral quartet calculation
#read in species tree
species_tree <- read.tree('concord.cf.tree')
#read in gene trees
gene_trees <- read.tree('Ophiostomatales_rRNA_PCGs_LOCI_nosp_nocontam.treefile')
genes <- c("SSU","ITS","LSU","Btub","EF1a","RPB1","RPB2")
for(i in 1:length(genes)){
  write.tree(keep.tip(species_tree,gene_trees[[i]]$tip.label), file = paste(genes[i],"species_tree.tre",sep = '_'))
}
