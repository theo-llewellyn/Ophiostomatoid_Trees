#!/bin/bash

#run ITSx to separate sequences into ITS, LSU and SSU
~/bin/ITSx_1.1.3/ITSx -i rRNA_operon.fasta --save_regions all --complement F --multi_thread T --cpu 8 --graphical F --preserve T --not_found F
