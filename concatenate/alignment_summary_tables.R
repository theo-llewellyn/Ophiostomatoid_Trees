library(tidyverse)

setwd("MAFFT/selected_representatives_woutgroups_contamremoved")

#list of genes
genes <- c('p1_SSU','p2_ITS','p3_LSU','p4_Btub','p5_EF1a','p6_RPB1','p7_RPB2')

#function to read in amas summary file, add a presence column and subset to just to taxon and presence
gene_function <- function(gene){
  table <- read_tsv(paste('Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp_nocontam_',gene,'-out.fas-seq-summary.txt', sep = ''))
  table[[gene]] <- rep(1,nrow(table))
  table <- table[,c(2,27)]
}
#apply to all genes
gene_tables <- lapply(genes, gene_function)

#join tables by taxon to make presence absence mat and then sum across columns
joined_table <- full_join(gene_tables[[1]], 
                          gene_tables[[2]], by = 'Taxon_name') %>% 
  full_join(., gene_tables[[3]]) %>%
  full_join(., gene_tables[[4]]) %>%
  full_join(., gene_tables[[5]]) %>%
  full_join(., gene_tables[[6]]) %>%
  full_join(., gene_tables[[7]]) %>%
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(across(where(is.numeric))))

#save file
write_tsv(joined_table, 'Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp_nocontam_summary_file.tsv')

