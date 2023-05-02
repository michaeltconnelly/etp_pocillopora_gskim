#!/bin/bash
#./bash/sereadconcat_wrapper.sh
#purpose: concatenate SE reads for spades assembly

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
cat ${prodir}/data/trimmed/${sample}_R1_SE_trimmed.fastq.gz ${prodir}/data/trimmed/${sample}_R2_SE_trimmed.fastq.gz > ${prodir}/data/trimmed/${sample}_SE_concat.fastq.gz
done