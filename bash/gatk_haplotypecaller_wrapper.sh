#!/bin/bash
#./bash/gatk_haplotypecaller_wrapper.sh
#purpose: call variants from processed bam files with GATK HaplotypeCaller

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

#loop to automate generation of scripts to call variants
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThC.q
#$ -l mres=64G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N gatk_HC_${sample}
#$ -o ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "module load bioinformatics/gatk" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "#" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
#
echo 'echo "This is the sample being processed:"' >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
echo "echo $sample" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
#   input command for GATK
echo "java -Xmx2g -jar /share/apps/bioinformatics/gatk/3.8.1.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-I ${prodir}/outputs/alignments/${sample}.sorted.md.rg.bam \
-o ${prodir}/outputs/alignments/${sample}.g.vcf.gz \
-R ${mcs}/sequences/pdam/pdam_genome.fasta \
-ERC GVCF \
-nct 16" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
#
echo "#" >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
# submit job
qsub ${prodir}/bash/jobs/gatk_haplotypecaller_${sample}.job
#
done