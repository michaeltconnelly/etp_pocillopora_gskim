#!/bin/bash
#./bash/oury_pipeline_align_wrapper.sh
#purpose: alignment to pocillopora UCE reference contigs

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

#loop to automate generation of scripts
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "#$ -j y
#$ -N ${sample}_oury_pipeline
#$ -o ${prodir}/bash/jobs/${sample}_oury_pipeline.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "module load bioinformatics/bwa" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "module load bioinformatics/samtools" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "module load bioinformatics/picard-tools" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "module load bioinformatics/gatk" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

#specify variables for program PATH
PICARD="/share/apps/bioinformatics/picard-tools/2.20.6/picard.jar"
GATK_HOME="/share/apps/bioinformatics/gatk/3.8.1.0"

#Example of clean reads’ filenames:
#${sample}_S1_R1_001_val_1.fq.gz	${sample}_S1_R2_001_val_2.fq.gz
#SA0002_S2_R1_001_val_1.fq.gz	SA0002_S2_R2_001_val_2.fq.gz

#Read mapping
#2.1°) Preparing the reference sequence - do this before starting to run wrapper jobs!
#Generate the BWA index
# bwa index ${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta
#Generate the fasta file index
# samtools faidx ${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta
#Generate the sequence dictionary
# java -Xmx4g -jar $PICARD CreateSequenceDictionary \
# REFERENCE=${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
# OUTPUT=${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.dict

echo "#2.2°) Mapping to the reference
#Decompressing fastq files" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "gzip -d ${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
echo "gzip -d ${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "#Mapping" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
# read group header
RG="@RG\tID:${sample}.XXXXXXXXX.1\tSM:${sample}\tPL:illumina\tLB:${sample}.S1\tPU:XXXXXXXXX.1"
#
echo "bwa mem \
-R '"$RG"' \
${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq \
${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq \
> ${prodir}/outputs/oury_pipeline/${sample}_aligned_reads.sam" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "#2.3°) Duplicate marking
#Sorting and converting to BAM" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "java -Xmx4g -jar $PICARD SortSam \
I=${prodir}/outputs/oury_pipeline/${sample}_aligned_reads.sam \
O=${prodir}/outputs/oury_pipeline/${sample}_sorted_reads.bam \
SORT_ORDER=coordinate" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "#Marking duplicates" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "java -Xmx4g -jar $PICARD MarkDuplicates \
I=${prodir}/outputs/oury_pipeline/${sample}_sorted_reads.bam \
O=${prodir}/outputs/oury_pipeline/${sample}_dedup_reads.bam \
M=${prodir}/outputs/oury_pipeline/${sample}_dup_metrics.txt \
ASO=coordinate" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "java -Xmx4g -jar $PICARD BuildBamIndex \
I=${prodir}/outputs/oury_pipeline/${sample}_dedup_reads.bam" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "#2.4°) Local realignment around indels" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "java -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R ${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
-I ${prodir}/outputs/oury_pipeline/${sample}_dedup_reads.bam \
-o ${prodir}/outputs/oury_pipeline/${sample}_target_intervals.list" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "java -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner \
-R ${mcs}/sequences/Pocillopora_UCEs_exons_2068_ref_sequences.fasta \
-I ${prodir}/outputs/oury_pipeline/${sample}_dedup_reads.bam \
-targetIntervals ${prodir}/outputs/oury_pipeline/${sample}_target_intervals.list \
-o  ${prodir}/outputs/oury_pipeline/${sample}_realigned_reads.bam" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

echo "#2.5°) Mapping QC" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo "samtools flagstat ${prodir}/outputs/oury_pipeline/${sample}_realigned_reads.bam" >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

# Echo finished 
echo 'echo "'${sample}' Finished!"' >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_oury_pipeline.job

# Submit job
qsub ${prodir}/bash/jobs/${sample}_oury_pipeline.job
#
done
