#!/bin/bash
#./bash/zip_wrapper.sh

set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples to be zipped:"
echo $samples

#loop to automate generation of scripts
for sample in $samples
do \
echo "gzip ${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq" >> ${prodir}/bash/jobs/${sample}_zip.job
echo "gzip ${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq" >> ${prodir}/bash/jobs/${sample}_zip.job

# Echo finished 
echo 'echo "'${sample}' Finished!"' >> ${prodir}/bash/jobs/${sample}_zip.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_zip.job

# Submit job
qsub ${prodir}/bash/jobs/${sample}_zip.job
#
done
