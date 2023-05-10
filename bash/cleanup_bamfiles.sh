#!/bin/bash
#./bash/cleanup_bamfiles.sh
# Cleanup script to manage hydra disk usage

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples with intermediate files to be removed:"
echo $samples

# remove intermediate bam files created by samtools and picard processing steps for ANGSD
for sample in $samples
do \
echo "removing $sample intermediate files"
rm ${prodir}/outputs/alignments/${sample}.sam
rm ${prodir}/outputs/alignments/${sample}.bam
rm ${prodir}/outputs/alignments/${sample}.sorted.bam
rm ${prodir}/outputs/alignments/${sample}.sorted.rg.bam
done
