#!/bin/bash

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"


for i in ${prodir}/outputs/mitofinder/*/*Final*/*contig.fasta ; do echo $(basename $i) ; grep -v ">" $i | wc -c ; done

for i in ${prodir}/outputs/mitofinder/*/*Final*/*contig.fasta ; do cat $i >> ${prodir}/outputs/mitofinder/all_complete_mtDNA_contigs.fasta ; done
