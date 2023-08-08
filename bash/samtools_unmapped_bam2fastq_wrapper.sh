#!/bin/bash
#./bash/samtools_unmapped_bam2fastq_wrapper.sh
#purpose: extract unmapped paired reads from .bam files with samtools

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples to be processed:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
echo "#$ -j y
#$ -N ${sample}_unmapped_bam2fastq
#$ -o ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
echo "module load bioinformatics/samtools" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job

echo 'echo "Starting samtools bam conversion to fastq"' >> $prodir/bash/jobs/${sample}_unmapped_bam2fastq.job

#   input command for samtools conversion
# -f 12 -F 256 outputs only reads where both reads of the pair are unmapped
echo "samtools bam2fq \
 -f 12 -F 256 \
 -1 ${prodir}/data/unmapped/${sample}_Unmapped_R1_PE.fastq \
 -2 ${prodir}/data/unmapped/${sample}_Unmapped_R2_PE.fastq \
 -@ 8 \
 ${prodir}/outputs/alignments/${sample}.bam" >> "${prodir}"/bash/jobs/${sample}_unmapped_bam2fastq.job

 echo "gzip ${prodir}/data/unmapped/${sample}_Unmapped_*.fastq" >> "${prodir}"/bash/jobs/${sample}_unmapped_bam2fastq.job

 echo 'echo '${sample}' successfully processed' >> "${prodir}"/bash/jobs/${sample}_unmapped_bam2fastq.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_unmapped_bam2fastq.job

 done