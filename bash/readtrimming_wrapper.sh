#!/bin/bash
#./bash/tagseq/readtrimming_wrapper.sh
#purpose: create wrapper scripts for read trimming using bbduk from BBtools suite

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/rsmas_samples.txt)

#lets me know which files are being processed
echo "These are the samples to be trimmed:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_trim.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=64G,h_data=4G,h_vmem=4G" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "#$ -j y
#$ -N ${sample}_trim
#$ -o ${prodir}/bash/jobs/${sample}_trim.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "module load bioinformatics/trimmomatic" >> ${prodir}/bash/jobs/${sample}_trim.job
echo 'module load java/1.8' >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_trim.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_trim.job

#   input command to trim raw reads
echo 'echo 'Trimming ${sample}'' >> "${prodir}"/bash/jobs/${sample}_trim.job
#
echo 'trimmomatic \
PE \
${prodir}/data/raw/${sample}_R1_001.fastq.gz \
${prodir}/data/raw/${sample}_R2_001.fastq.gz \
${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz \
${prodir}/data/trimmed/${sample}_R1_SE_trimmed.fastq.gz \
${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz \
${prodir}/data/trimmed/${sample}_R2_SE_trimmed.fastq.gz \
ILLUMINACLIP:${prodir}/data/itru_adapters_trimmomatic.fa:2:20:5:3:keepBothReads \
CROP:150 \
HEADCROP:3 \
SLIDINGWINDOW:5:15 \
MINLEN:30' >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo 'echo '${sample}' successfully trimmed' >> "${prodir}"/bash/jobs/${sample}_trim.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_trim.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_trim.job
#
done
