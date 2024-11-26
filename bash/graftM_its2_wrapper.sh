#!/bin/bash
#./bash/graftM_its2_wrapper.sh
#purpose: perform graftM assignment of ITS2 sequences from non-coral reads

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples to be assembled:"
echo $samples

# loop to automate generation of scripts to direct sequence file trimming
# edit output job script location within your project root directory
JOBFILE="${prodir}/bash/jobs/${sample}_graftM_its2.job"

echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=16G,h_data=4G,h_vmem=4G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_graftM_its2
#$ -o ${prodir}/bash/jobs/${sample}_graftM_its2.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE

echo "# ----------------Modules------------------------- #" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

for sample in $samples
do \
# insert graftM command
echo "graftM graft --forward data/unmapped/${sample}_Unmapped_R1_PE.fastq.gz \
--reverse data/unmapped/${sample}_Unmapped_R2_PE.fastq.gz \
--graftm_package ${prodir}/data/seqs/ITS2_graftm_final.gpkg \
--input_sequence_type nucleotide \
--output_directory outputs/graftM/${sample}" >> $JOBFILE

echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
done

