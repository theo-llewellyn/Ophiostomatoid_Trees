#!/bin/bash
#subset TBAS 6-loci fasta files data for Magnaporthales 
for file in *.fas;
do seqtk subseq $file outgroups.txt > outgroups_${file};
done

#align outgroups to existing alignments
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_SSU.fas --thread 8 selected_representatives/Ophiostomatales_selected_sequences_nosp_p1_SSU-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p1_SSU-out.fas
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_ITS.fas --thread 8 selected_representatives/Ophiostomatales_selected_sequences_nosp_p2_ITS-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p2_ITS-out.fas
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_LSU.fas --thread 1 selected_representatives/Ophiostomatales_selected_sequences_nosp_p3_LSU-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p3_LSU-out.fas
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_EF1a.fas --thread 8 selected_representatives/Ophiostomatales_selected_sequences_nosp_p5_EF1a-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p5_EF1a-out.fas
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_RPB1.fas --thread 8 selected_representatives/Ophiostomatales_selected_sequences_nosp_p6_RPB1-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p6_RPB1-out.fas
mafft --auto --addfragments ../outgroup_nexusFiles/outgroups_RPB2.fas --thread 8 selected_representatives/Ophiostomatales_selected_sequences_nosp_p7_RPB2-out.fas > selected_representatives_woutgroups/Ophiostomatales_selected_sequences_nosp_p7_RPB2-out.fas

# concatenate
~/bin/AMAS/amas/AMAS.py concat -i Ophiostomatales_selected_sequences_p*-out.fas -f fasta -d dna -t Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups.fasta -p Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.txt
~/bin/AMAS/amas/AMAS.py concat -i Ophiostomatales_selected_sequences_nosp*-out.fas -f fasta -d dna -t Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp.fasta -p Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.nosp.txt

#Update summary tables showing which sequences have which genes
gsed -i 's/p.*_Ophiostomatales_selected_sequences_nosp_p*._//;s/-out//;s/DNA, //' Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.nosp.txt
~/bin/AMAS/amas/AMAS.py split -i Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_nosp_nocontam.fasta -f fasta -d dna -l Ophiostomatales_rRNA_PCGs_concat_msa_woutgroups_partition.nosp.txt -j
~/bin/AMAS/amas/AMAS.py summary -f fasta -d dna -i Ophiostomatales_rRNA_PCGs_concat_msa_*-out.fas -s
