#!/bin/bash

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# print filename, length of mitogenome
for i in ${prodir}/outputs/mitofinder/*/*Final*/*contig.fasta ; do echo $(basename $i) ; grep -v ">" $i | wc -c ; done

# combine complete contigs into a single fasta file
for i in ${prodir}/outputs/mitofinder/*/*Final*/*contig.fasta ; do cat $i >> ${prodir}/outputs/mitofinder/all_complete_mtDNA_contigs.fasta ; done
