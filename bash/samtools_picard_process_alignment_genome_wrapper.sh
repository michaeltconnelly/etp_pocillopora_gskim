#!/bin/bash
#./bash/samtools_picard_process_alignment_wrapper.sh
#purpose: sort, index, and convert alignments to .bam files with samtools, add read groups and mark duplicates with picard

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

# specify reference genome
genome=$2

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo "#$ -j y
#$ -N ${sample}_samtools_picard
#$ -o ${prodir}/bash/jobs/${sample}_samtools_picard.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo "module load bioinformatics/samtools" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo "module load bioinformatics/picard-tools" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_samtools_picard.job

echo 'echo "Starting samtools bam conversion, sort and index steps"' >> $prodir/bash/jobs/${sample}_samtools_picard.job
#   input command for samtools conversion
echo "samtools view -b ${prodir}/outputs/alignments/${genome}/${sample}.sam \
-o ${prodir}/outputs/alignments/${genome}/${sample}.bam -@ 8"  >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
#
echo "samtools sort \
${prodir}/outputs/alignments/${genome}/${sample}.bam -@ 8 \
-o ${prodir}/outputs/alignments/${genome}/${sample}.sorted.bam" >> $prodir/bash/jobs/${sample}_samtools_picard.job
echo "#" >> $prodir/bash/jobs/${sample}_samtools_picard.job
#
echo "samtools index -b \
${prodir}/outputs/alignments/${genome}/${sample}.sorted.bam" >> $prodir/bash/jobs/${sample}_samtools_picard.job
echo "#" >> $prodir/bash/jobs/${sample}_samtools_picard.job

#   input command for picard add read groups 
echo 'echo "Starting PicardTools and GATK processing steps"' >> $prodir/bash/jobs/${sample}_samtools_picard.job
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
AddOrReplaceReadGroups \
INPUT=${prodir}/outputs/alignments/${genome}/${sample}.sorted.bam \
OUTPUT=${prodir}/outputs/alignments/${genome}/${sample}.sorted.rg.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=unit1 \
RGSM=${sample}" >> $prodir/bash/jobs/${sample}_samtools_picard.job
echo "#" >> $prodir/bash/jobs/${sample}_samtools_picard.job
#
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
MarkDuplicates \
INPUT=${prodir}/outputs/alignments/${genome}/${sample}.sorted.rg.bam \
OUTPUT=${prodir}/outputs/alignments/${genome}/${sample}.sorted.md.rg.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=${prodir}/outputs/alignments/${genome}/${sample}_marked_dup_metrics.txt" >> $prodir/bash/jobs/${sample}_samtools_picard.job
echo "#" >> $prodir/bash/jobs/${sample}_samtools_picard.job
#
echo 'echo '${sample}' successfully processed' >> "${prodir}"/bash/jobs/${sample}_samtools_picard.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_samtools_picard.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_samtools_picard.job
#
sleep 0.2
done