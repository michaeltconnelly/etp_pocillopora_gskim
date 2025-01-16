#!/bin/bash
#./bash/samtools_symbiontCD_process_alignment_wrapper.sh
#purpose: sort, index, and convert symbiont_alignments to .bam files with samtools, add read groups and mark duplicates with 

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
JOBFILE="${prodir}/bash/jobs/${sample}_symCD_samtools.job"
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
#$ -N ${sample}_symCD_samtools
#$ -o ${prodir}/bash/jobs/${sample}_symCD_samtools.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/samtools" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

echo 'echo "Starting samtools bam conversion, sort and index steps"' >> $JOBFILE

#   input command for samtools conversion to bam, produce separate bam files for Cladocopium and Durusdinium alignments
# Cladocopium (chr11)
echo "samtools view -b ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sam \
-o ${prodir}/outputs/symbiont_alignments/${sample}_symC.bam -@ 8 chr11"  >> $JOBFILE
# Durusdinium (chr12) 
echo "samtools view -b ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sam \
-o ${prodir}/outputs/symbiont_alignments/${sample}_symD.bam -@ 8 chr12"  >> $JOBFILE

#   input command for samtools sort and index 
# Cladocopium
echo "samtools sort \
${prodir}/outputs/symbiont_alignments/${sample}_symC.bam -@ 8 \
> ${prodir}/outputs/symbiont_alignments/${sample}_symC.sorted.bam" >> $JOBFILE
#
echo "samtools index -b \
${prodir}/outputs/symbiont_alignments/${sample}_symC.sorted.bam" >> $JOBFILE
echo "#" >> $JOBFILE
# Durusdinium
echo "samtools sort \
${prodir}/outputs/symbiont_alignments/${sample}_symD.bam -@ 8 \
> ${prodir}/outputs/symbiont_alignments/${sample}_symD.sorted.bam" >> $JOBFILE
#
echo "samtools index -b \
${prodir}/outputs/symbiont_alignments/${sample}_symD.sorted.bam" >> $JOBFILE
echo "#" >> $JOBFILE

echo 'echo '${sample}' successfully processed' >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
sleep 0.2
done