#!/bin/bash
#./bash/bwa_align_filter_contaminants_wrapper.sh
#purpose: alignment to pocillopora reference genome

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
# create job file
JOBFILE="${prodir}/bash/jobs/${sample}_bwa_decontam_align.job"
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
#$ -N ${sample}_bwa_align
#$ -o ${prodir}/bash/jobs/${sample}_bwa_decontam_align.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- # 
module load bioinformatics/bwa
module load bioinformatics/samtools" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

  #   input bwa mem alignment command
echo "bwa mem ${prodir}/data/seqs/sym_mtDNA/sym_mtDNA_seqs.fasta \
${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz \
${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz \
> ${prodir}/outputs/contaminant_alignments/${sample}.sam" >> $JOBFILE

  #   input command for samtools conversion
echo "samtools view -b ${prodir}/outputs/contaminant_alignments/${sample}.sam \
-o ${prodir}/outputs/contaminant_alignments/${sample}.bam -@ 8" >> $$JOBFILE

echo "#" >> $$JOBFILE

  #   input command to filter out aligned contaminant reads
echo 'echo "Starting samtools bam conversion to fastq"' >> $JOBFILE

  #   input command for samtools conversion
  # -f 12 -F 256 outputs only reads where both reads of the pair are unmapped
echo "samtools bam2fq \
 -f 12 -F 256 \
 -1 ${prodir}/data/decontaminated/${sample}_decontam_R1_PE.fastq \
 -2 ${prodir}/data/decontaminated/${sample}_decontam_R2_PE.fastq \
 -@ 8 \
 ${prodir}/outputs/contaminant_alignments/${sample}.bam" >> $JOBFILE

echo "#" >> $$JOBFILE

echo "gzip ${prodir}/data/decontaminated/${sample}_decontam_*.fastq" >> $JOBFILE
#
echo 'echo '${sample}' successfully processed' >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
done