#remove spaces and dots in msa to match the names changed by AMAS
sed "s/ /_/g;s/\./_/g;s/\'//g" Ophiostomatales_rRNA_PCGs_concat_msa.fasta > Ophiostomatales_rRNA_PCGs_concat_renamed.msa.fasta

#select best of each species, using first gene coverage and then percentage gaps missing to break ties
python sequence_selecter_new.py

cp Ophiostomatales_rRNA_PGCs_concat_msa_partition_ITSas1.txt Ophiostomatales_selected_sequences.tsv Ophiostomatales_selected_sequences.fasta > selected_representatives
cd selected_representatives

#separate into each gene
~/bin/AMAS/amas/AMAS.py split -i Ophiostomatales_selected_sequences.fasta -f fasta -d dna -l Ophiostomatales_rRNA_PGCs_concat_msa_partition_ITSas1.txt -j 

#remove gap columns and reconcatenate
~/bin/AMAS/amas/AMAS.py concat -i *-out.fas -f fasta -d dna -t Ophiostomatales_selected_sequences.fasta -p Ophiostomatales_rRNA_PCGs_concat_msa_partition.txt

#make copy with only sequences with full binomial, i.e. no sp. will still keep things from collections e.g. Ceratocystis_CBSXXX
grep '_sp' Ophiostomatales_selected_sequences.fasta | grep -v '_sp[[:alnum:]]*_' | sed 's/>//g' > sp_to_remove.txt
seqkit grep -v -f sp_to_remove.txt Ophiostomatales_selected_sequences.fasta > Ophiostomatales_selected_sequences_nosp.fasta

#split into genes
~/bin/AMAS/amas/AMAS.py split -i Ophiostomatales_selected_sequences_nosp.fasta -f fasta -d dna -l Ophiostomatales_rRNA_PCGs_concat_msa_partition.txt -j 
