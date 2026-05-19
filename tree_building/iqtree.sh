#!/bin/bash
iqtree \
 -alrt 1000 -abayes -lbp 1000 -bnni -B 1000 \
 -s ../01_ALIGNMENT/Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp_nocontam.fasta \
 -nt AUTO --threads-max 32 \
 -cptime 120 \
 --prefix Ophiostomatales_rRNA_PCGs_nosp_nocontam \
 -p ../01_ALIGNMENT/Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.nosp.txt
