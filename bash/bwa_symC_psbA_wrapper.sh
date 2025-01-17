#!/bin/bash
#./bash/bwa_symC_psbA_wrapper.sh
#purpose: alignment to Cladocopium psbA reference and consensus sequence extraction

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
JOBFILE="${prodir}/bash/jobs/${sample}_bwa_symC_psbA.job"
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_bwa_symC_psbA
#$ -o ${prodir}/bash/jobs/${sample}_bwa_symC_psbA.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/bwa" >> $JOBFILE
echo "module load bioinformatics/samtools" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

#   pipeline modified from: https://github.com/hisatakeishida/Symb-SHIN/blob/main/C_Read-based_alignment.md

#   input bwa mem alignment command
echo "bwa mem ${prodir}/data/seqs/cladocopium_psbA.fasta \
${prodir}/data/unmapped/${sample}_Unmapped_R1_PE.fastq.gz \
${prodir}/data/unmapped/${sample}_Unmapped_R2_PE.fastq.gz \
> ${prodir}/outputs/symC_psbA/${sample}_C_psbA.sam" >> $JOBFILE

#   input convert to bam command
echo "samtools view -b -@ 8 -F 4 \
${prodir}/outputs/symC_psbA/${sample}_C_psbA.sam \
-o ${prodir}/outputs/symC_psbA/${sample}_C_psbA.bam" >> $JOBFILE

#   keep only full-length reads aligned with more than 90% of identity
echo "perl ${prodir}/bash/BAMFiltration.pl -in ${prodir}/outputs/symC_psbA/${sample}_C_psbA.bam \
 -out ${prodir}/outputs/symC_psbA/${sample}_C_psbA.filtered100-90.bam -minPCaligned 100 -minPCidentity 90" >> $JOBFILE

echo 'echo '${sample}' successfully processed' >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
done

# proceed with identifying the reference with highest number of reads mapped, generating text mpileup output and building consensus sequences
# ./bash/symC_psbA_consensus_seqs.sh