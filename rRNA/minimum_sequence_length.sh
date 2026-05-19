#!/bin/bash

#min 100bp for SSU and LSU
awk 'BEGIN {RS=">"; ORS=""} 
     NR>1 {
       n=0; seq=""; 
       split($0, lines, "\n"); 
       header=lines[1]; 
       for(i=2;i<=length(lines);i++) seq=seq lines[i]; 
       gsub(/[ \t\r\n]/,"",seq); 
       if(length(seq) >= 100) print ">"header"\n"seq"\n"
     }' ITSx_out.SSU.fasta > ITSx_out.SSU.min100bp.fasta
     
awk 'BEGIN {RS=">"; ORS=""} 
     NR>1 {
       n=0; seq=""; 
       split($0, lines, "\n"); 
       header=lines[1]; 
       for(i=2;i<=length(lines);i++) seq=seq lines[i]; 
       gsub(/[ \t\r\n]/,"",seq); 
       if(length(seq) >= 100) print ">"header"\n"seq"\n"
     }' ITSx_out.LSU.fasta > ITSx_out.LSU.min100bp.fasta
