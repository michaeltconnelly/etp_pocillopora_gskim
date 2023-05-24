#!/bin/bash
# purpose: extract the CDS for the 13 protein-coding genes and 2 rRNA subunits annotated in Mitofinder results

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# extract CDS from the NT fasta files
mtgenes="ATP8 COX1 rrnL ND1 CYTB ND2 ND6 ATP6 ND4 rrnS COX3 COX2 ND4L ND3 ND5"

# prep directories and output files --> add flow control if statements
for gene in $mtgenes; do mkdir ${prodir}/outputs/mitofinder/cds_alignments/${gene} 
touch ${prodir}/outputs/mitofinder/cds_alignments/${gene}/${gene}_all.fasta; done

# nested for-loops to extract CDS from each sample into single file for each gene
# complete assemblies
for gene in $mtgenes; do for file in ${prodir}/outputs/mitofinder/final_results/*/*contig_genes_NT.fasta; do grep -A 1 ">*\@${gene}" $file >> ${prodir}/outputs/mitofinder/cds_alignments/${gene}/${gene}_all.fasta; done
done
# incomplete assemblies
for gene in $mtgenes; do for file in ${prodir}/outputs/mitofinder/final_results/*/*final_genes_NT.fasta; do grep -A 1 ">*\@${gene}" $file >> ${prodir}/outputs/mitofinder/cds_alignments/${gene}/${gene}_all.fasta; done
done

# print message
echo "MitoFinder CDS extraction and concatenation complete"
for gene in $mtgenes; do for i in ${prodir}/outputs/mitofinder/cds_alignments/${gene}/${gene}_all.fasta ; do lines=$(grep -v ">" $i | wc -l) ; echo "There are $lines seqs in $(basename $i)" ; done ; done
#