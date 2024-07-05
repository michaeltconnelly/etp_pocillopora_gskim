#!/bin/bash
#./bash/gatk_combinegvcfs_wrapper.sh
#purpose: combine per-sample gvcf files from GATK HaplotypeCaller into a cohort gvcf file

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

#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q lThC.q
#$ -l mres=64G,h_data=8G,h_vmem=8G
#$ -cwd
#$ -j y
#$ -N gatk_combine_gvcfs
#$ -o gatk_combine_gvcfs.log
#$ -m bea
#$ -M connellym@si.edu
#"  >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo "module load bioinformatics/gatk" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/gatk_combinegvcfs.job

# Input command for GATK CombineGVCFs
echo "java -jar -jar /share/apps/bioinformatics/gatk/4.5.0.0/gatk-package-4.5.0.0-local.jar \
-T CombineGVCFs \
-R ${mcs}/sequences/p_grandis_GCA_964027065.2/GCA_964027065.2_jaPocGran1.hap1.2_genomic.fna \\" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
# for loop to input which gvcfs to combine
for sample in $samples
do \
echo "--variant ${prodir}/outputs/alignments/gvcfs/${sample}.g.vcf.gz \\" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
done
# 
echo "--out ${prodir}/outputs/${set}_g.vcf.gz" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
#
echo "#" >> ${prodir}/bash/jobs/gatk_combinegvcfs.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/gatk_combinegvcfs.job

# submit job
qsub ${prodir}/bash/jobs/gatk_combinegvcfs.job