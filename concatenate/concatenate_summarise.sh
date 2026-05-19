#concatenate alignments
~/bin/AMAS/amas/AMAS.py concat -i \
 SSU.renamed.msa.fasta \
 ITS1.renamed.msa.fasta \
 5.8S.renamed.msa.fasta \
 ITS2.renamed.msa.fasta \
 LSU.renamed.msa.fasta \
 Btub.renamed.msa.fasta \
 EF1a.renamed.msa.fasta \
 RPB1.renamed.msa.fasta \
 RPB2.renamed.msa.fasta \
 -f fasta -d dna \
 -t Ophiostomatales_rRNA_PCGs_concat_msa.fasta \
 -p Ophiostomatales_rRNA_PGCs_concat_msa_partition.txt


#split the rRNA alignment with ITS as a single partition and get summary
~/bin/AMAS/amas/AMAS.py split -i Ophiostomatales_rRNA_PCGs_concat_msa.fasta -f fasta -d dna -l Ophiostomatales_rRNA_PGCs_concat_msa_partition_ITSas1.txt -j 

~/bin/AMAS/amas/AMAS.py summary -f fasta -d dna -i Ophiostomatales_rRNA_PCGs_concat_msa_p*-out.fas -s
