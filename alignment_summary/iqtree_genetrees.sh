#!/bin/bash
#make gene trees
iqtree -s ../01_ALIGNMENT/Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp_nocontam.fasta \
 -S ../01_ALIGNMENT/Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.nosp.txt \
 --prefix Ophiostomatales_rRNA_PCGs_LOCI_nosp_nocontam \
 -T AUTO \
 -B 1000 \
 -bnni \
 --threads-max 32 \
 --redo

#calculate gene and site concordance
iqtree -t Ophiostomatales_rRNA_PCGs_nosp_nocontam.treefile \
 --gcf Ophiostomatales_rRNA_PCGs_LOCI_nosp_nocontam.treefile \
 -p ../01_ALIGNMENT/forGCF \
 --scf 100 --prefix concord
