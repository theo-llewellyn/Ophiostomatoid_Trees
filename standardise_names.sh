#!/bin/bash
#make a table from taxa headers with the sequence name and sample information
loci=(RPB1 RPB2 EF1a Btub SSU ITS1 5.8S ITS2 LSU)
for locus in $loci
do
awk 'BEGIN{OFS="\t"}
  /^>/ {
    header = substr($0,2)                     # remove leading ">"
    species = $2 " " $3
    rest = substr(header, index(header,$2))   # skip accession

    # extract strain/isolate/culture IDs or codes
    cmd = "echo \"" rest "\" | grep -oE \"" \
          "strain[: ]+[A-Za-z0-9 _.-]+|" \
          "isolate[: ]+[A-Za-z0-9 _.-]+|" \
          "culture[: ]+[A-Za-z0-9 _:<>._-]+|" \
          "\\([^)]+\\)|" \
          "[A-Z]{1,5}(-[A-Z0-9]+)+|" \
          "[A-Z]{2,}[ ]+[0-9]{1,5}|" \
          "(KM|CWM|CBS)[ ]*[A-Za-z0-9 _-]*\""

    cmd | getline matches
    close(cmd)

    gsub(/^ +| +$/,"",matches)   # trim whitespace
    gsub(/^(strain[: ]+|isolate[: ]+|culture[: ]+)/,"",matches)  # remove prefixes

    if(matches == "") matches = ""   # ensures empty second column if no match

    print header, species, matches
  }' ${locus}_msa.fasta > ${locus}_sequence_strain_table.tsv
done

#remove any punctuation and spaces in voucher codes using the follow regular expression
[[:punct:] ]
#remove / from species name column and make seds in excel using the following excel command
="s/"&A1&"/"&B1&" "&C1&"/g"

#replace names in fasta files to species and then code
for locus in $loci
do
sed -f ${locus}_seds.txt ${locus}_msa.fasta > ${locus}.renamed.msa.fasta
done

#standardise sequences labelled as sp.
for file in *.renamed.msa.fasta;
do 
gsed -i -f final_seds.txt $file 
gsed -i 's/sp. /sp/g;s/sp./sp/g' $file;
done
