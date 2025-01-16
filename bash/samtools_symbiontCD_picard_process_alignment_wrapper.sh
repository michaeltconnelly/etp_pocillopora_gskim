#!/bin/bash
#./bash/samtools_symbiontCD_picard_process_alignment_wrapper.sh
#purpose: sort, index, and convert symbiont_alignments to .bam files with samtools, add read groups and mark duplicates with picard

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
JOBFILE="${prodir}/bash/jobs/${sample}_symCD_samtools_picard.job"
touch $JOBFILE
#
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_symCD_samtools_picard
#$ -o ${prodir}/bash/jobs/${sample}_symCD_samtools_picard.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/samtools" >> $JOBFILE
echo "module load bioinformatics/picard-tools" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

echo 'echo "Starting samtools bam conversion, sort and index steps"' >> $JOBFILE
#   input command for samtools conversion
echo "samtools view -b ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sam \
-o ${prodir}/outputs/symbiont_alignments/${sample}_symCD.bam -@ 8"  >> $JOBFILE
#
echo "samtools sort \
${prodir}/outputs/symbiont_alignments/${sample}_symCD.bam -@ 8 \
> ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.bam" >> $JOBFILE
echo "#" >> $JOBFILE
#
echo "samtools index -b \
${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.bam" >> $JOBFILE
echo "#" >> $JOBFILE

#   input command for picard add read groups 
echo 'echo "Starting PicardTools and GATK processing steps"' >> $JOBFILE
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
AddOrReplaceReadGroups \
INPUT=${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.bam \
OUTPUT=${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.rg.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=unit1 \
RGSM=${sample}" >> $JOBFILE
echo "#" >> $JOBFILE
#
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
MarkDuplicates \
INPUT=${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.rg.bam \
OUTPUT=${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.md.rg.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=${prodir}/outputs/symbiont_alignments/${sample}_symCD_marked_dup_metrics.txt" >> $JOBFILE
echo "#" >> $JOBFILE
#

#   input commands to produce separate bam files for Cladocopium and Durusdinium alignments
# Cladocopium (chr11)
echo "samtools view -b ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.md.rg.bam \
-o ${prodir}/outputs/symbiont_alignments/${sample}_symC.sorted.md.rg.bam -@ 8 chr11"  >> $JOBFILE
# Durusdinium (chr12)
echo "samtools view -b ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sorted.md.rg.bam \
-o ${prodir}/outputs/symbiont_alignments/${sample}_symD.sorted.md.rg.bam -@ 8 chr12"  >> $JOBFILE

echo 'echo '${sample}' successfully processed' >> "${prodir}"/bash/jobs/${sample}_samtools_picard.job
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
sleep 0.2
done