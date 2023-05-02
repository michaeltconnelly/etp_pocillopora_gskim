#!/bin/bash
#./bash/spades_wrapper.sh
#purpose: spades assembly

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples to be trimmed:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_spades.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_spades.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mmres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_spades.job
echo "#$ -j y
#$ -N ${sample}_spades
#$ -o ${prodir}/bash/jobs/${sample}_spades.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_spades.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_spades.job
echo "module load bioinformatics/spades" >> ${prodir}/bash/jobs/${sample}_spades.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_spades.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_spades.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_spades.job

#   input command for spades assembly 
echo "spades.py \
-o ${prodir}/outputs/spades/${sample} \
--pe1-1 ${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz \
--pe1-2 ${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz \
--pe1-s ${prodir}/data/trimmed/${sample}_SE_concat.fastq.gz"  >> ${prodir}/bash/jobs/${sample}_spades.job
#
echo 'echo '${sample}' successfully assembled' >> "${prodir}"/bash/jobs/${sample}_spades.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_spades.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_spades.job
#
done