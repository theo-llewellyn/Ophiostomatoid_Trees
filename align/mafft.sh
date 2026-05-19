#!/bin/bash

#align
for i in *.fasta
do PREFIX=${i%%.fasta*}
echo ${PREFIX}
mafft --auto $i > ../MAFFT/${PREFIX}_msa.fasta
done
