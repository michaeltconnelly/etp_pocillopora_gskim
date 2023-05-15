#!/bin/bash
#./bash/bwa_symbiont_wrapper.sh
#purpose: alignment to concatenated symbiont reference genome

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
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
echo "#$ -j y
#$ -N ${sample}_bwa_align_symbiont
#$ -o ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
echo "module load bioinformatics/bwa" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job

  #   input bwa mem alignment command
echo "bwa mem ${mcs}/sequences/symbio/symABCDgenome.fasta \
${prodir}/data/unmapped/${sample}_Unmapped_R1_PE.fastq.gz \
${prodir}/data/unmapped/${sample}_Unmapped_R2_PE.fastq.gz \
> ${prodir}/outputs/alignments/${sample}_symABCD.sam" >> "${prodir}"/bash/jobs/${sample}_bwa_align_symbiont.job

#
echo 'echo '${sample}' successfully aligned' >> "${prodir}"/bash/jobs/${sample}_bwa_align_symbiont.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_bwa_align_symbiont.job
#
done