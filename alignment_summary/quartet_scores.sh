# List of genes
genes=("SSU" "ITS" "LSU" "Btub" "EF1a" "RPB1" "RPB2")

#save gene trees as separate files e.g.
head -n 7 Ophiostomatales_rRNA_PCGs_LOCI_nosp_nocontam.treefile | tail -n 1 > RPB2.treefile

# Output summary table
echo -e "Gene\tFinalNormalizedQuartetScore" > quartet_scores.tsv

# Loop through each gene
for gene in "${genes[@]}"; do
  echo "Running ASTRAL for $gene..."

  # Run ASTRAL
  java -jar ~/bin/Astral/astral.5.7.8.jar \
    -q IQTree/${gene}_species_tree.tre \
    -i IQTree/${gene}.treefile \
    -t 1 > ${gene}_output.txt 2>&1

  # Extract the "Final normalized quartet score" value
  score=$(grep "normalized quartet" ${gene}_output.txt | tail -n 1 | cut -f 6 -d ' ')

  # Write to table
  echo -e "${gene}\t${score}" >> quartet_scores.tsv
  rm ${gene}_output.txt
done

#all genes combined
java -jar ~/bin/Astral/astral.5.7.8.jar \
    -q IQTree/Ophiostomatales_rRNA_PCGs_nosp_nocontam.treefile \
    -i IQTree/Ophiostomatales_rRNA_PCGs_LOCI_nosp_nocontam.treefile \
    -t 1 2> Astral_output_all.txt 2>&1
score=$(grep "normalized quartet" Astral_output_all.txt | tail -n 1 | cut -f 6 -d ' ')

echo -e "7loci\t${score}" >> quartet_scores.tsv
