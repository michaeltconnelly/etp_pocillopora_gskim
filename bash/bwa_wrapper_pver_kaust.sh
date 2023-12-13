#!/bin/bash
#./bash/bwa_align_wrapper.sh
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
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
echo "#$ -j y
#$ -N ${sample}_pver-kaust_bwa_align
#$ -o ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
echo "module load bioinformatics/bwa" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job

  #   input bwa mem alignment command
echo "bwa mem ${mcs}/sequences/pver_kaust/GCA_014529365.1_Pver_genome_assembly_v1.0_genomic.fna \
${prodir}/data/trimmed/${sample}_R1_PE_trimmed.fastq.gz \
${prodir}/data/trimmed/${sample}_R2_PE_trimmed.fastq.gz \
> ${prodir}/outputs/alignments/pver_kaust/${sample}.sam" >> "${prodir}"/bash/jobs/${sample}_pver-kaust_bwa_align.job

#
echo 'echo '${sample}' successfully aligned' >> "${prodir}"/bash/jobs/${sample}_pver-kaust_bwa_align.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_pver-kaust_bwa_align.job
#
done