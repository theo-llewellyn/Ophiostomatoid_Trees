# A robust phylogenetic framework for ophiostomatoid fungi
This repository contains the code associated to the paper:
<br/>

A robust phylogenetic framework for ophiostomatoid fungi: orders Ophiostomatales and Microascales (Sordariomycetes, Ascomycota)
<br/>

**Authors**:

Theo Llewellyn<sup>1,2,</sup>*, Alfried Vogler<sup>1,2</sup>
<br/>

**Affilitions**<br/>
1. Leverhulme Centre for the Holobiont, Department of Life Sciences, Imperial College London, Silwood Park Campus, Ascot, Berkshire, SL5 7PY, UK
2. Department of Life Sciences, Natural History Museum, Cromwell Road, London, SW6 7BD UK

*Corresponding author: t.llewellyn19@imperial.ac.uk

## Data Records

The multiple sequence alignments and phylogenetic trees generated and analysed during the current study are deposited in the figshare repository: https://doi.org/10.6084/m9.figshare.30580976.v1 Figshare contains:
1. 
2.


## Analysis scripts
All bioinformatic code used to download and process sequences, build trees, and calculate summary statistics

### 1. Sequence Download
The search strings used to download sequences from NCBI can be found here: `./search_strings.txt`


### 2. rRNA operon split
`cd rRNA`
1. `./ITSx.sh` identifies and splits sequences into ITS, LSU and SSU
2. `./minimum_sequence_length.sh` sets 100bp minimum threshold for SSU and LSU

### 3. Alignment
`cd align`
1. `./mafft.sh`

### 4. Standardise fasta headers and concatenate alignments
`cd concatenate`
1. `./standardise_names.sh` standardise and rename taxa header to SPECIES_CULTURECODE
2. `./concatenate_summarise.sh` concatenate sequences and summarise in tables
3. `Rscript alignment_summary_tables.R` make table summarising which genes are present for which taxa
4. `./sequence_selecter.sh` select best sequence based on gene coverage and length. requires sequence_selecter_new.py python script.

### 5. Add outgroups
`cd outgroups`
1. `./outgroups.sh` adds outgroups to existing alignments, reconcatenates and updates summary tables

### 6. Tree building
`cd tree_building`
1. `./iqtree.sh` builds concatenated species tree and runs bootstrap analysis

### 7. Calculate alignment and gene summary statistics
1. `./iqtree_genetrees.sh` builds gene trees and calculate site and gene concordance compared to species tree
3. `./quartet_scores.sh` calculates quartet scores and outputs results in table
