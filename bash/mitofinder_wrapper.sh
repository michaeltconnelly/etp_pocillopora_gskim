#!/bin/bash
#./bash/mitofinder_wrapper.sh
#purpose: recover mitochondrial genome

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples for mitofinder:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_mitofinder.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
echo "#$ -j y
#$ -N ${sample}_mitofinder
#$ -o ${prodir}/bash/jobs/${sample}_mitofinder.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
echo "module load bioinformatics/mitofinder" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_mitofinder.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_mitofinder.job

#   input mitofinder command
# mitosequence.gb consists of Genbank reference database, 5749 mitogenomes 05/02/2023
# pocilloporidae.gb consists of Pocillopora (2), Madracis (1), Stylophora (1), Seriatopora (2) mitogenomes
# pocillopora.gb consists of Pocillopora-only reference mitogenomes
echo "mitofinder \
-j ${sample} \
-o 4 \
-r ${mcs}/sequences/pocillopora.gb \
-1 ${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz  \
-2 ${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz  \
--new-genes" >> "${prodir}"/bash/jobs/${sample}_mitofinder.job

#
echo 'echo '${sample}' successfully processed' >> "${prodir}"/bash/jobs/${sample}_mitofinder.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_mitofinder.job
# submit job
qsub -wd ${prodir}/outputs/mitofinder ${prodir}/bash/jobs/${sample}_mitofinder.job
#
done