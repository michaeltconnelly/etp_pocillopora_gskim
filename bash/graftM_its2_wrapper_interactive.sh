#!/bin/bash
#./bash/graftM_its2_wrapper_interactive.sh
#purpose: perform graftM assignment of ITS2 sequences from non-coral reads

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

for sample in $samples
do \

echo $sample

# graftM command
~/programs/graftM/bin/graftM graft --forward data/unmapped/${sample}_Unmapped_R1_PE.fastq.gz \
--reverse data/unmapped/${sample}_Unmapped_R2_PE.fastq.gz \
--graftm_package ${prodir}/data/seqs/ITS2_graftm_final.gpkg \
--input_sequence_type nucleotide \
--output_directory outputs/graftM/${sample}
done

