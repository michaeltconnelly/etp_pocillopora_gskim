#!/bin/bash
#./bash/bwa_symbiontCD_wrapper.sh
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
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
echo "#$ -j y
#$ -N ${sample}_bwa_align_symbiont_CD
#$ -o ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
echo "module load bioinformatics/bwa" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job

#   input bwa mem alignment command
# NOTE: use symCDgenome.fasta for downstream use with zooxtype.pl (all scaffolds collapsed to single header)
# NOTE: use sym_catCD_genome.fasta for original genomes (all scaffolds retained)
echo "bwa mem ${mcs}/sequences/symbio/symCDgenome.fasta \
${prodir}/data/unmapped/${sample}_Unmapped_R1_PE.fastq.gz \
${prodir}/data/unmapped/${sample}_Unmapped_R2_PE.fastq.gz \
> ${prodir}/outputs/symbiont_alignments/${sample}_symCD.sam" >> "${prodir}"/bash/jobs/${sample}_bwa_align_symbiont_CD.job

#
echo 'echo '${sample}' successfully aligned' >> "${prodir}"/bash/jobs/${sample}_bwa_align_symbiont_CD.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_bwa_align_symbiont_CD.job
#
done