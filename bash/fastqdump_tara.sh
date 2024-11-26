#!/bin/bash

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
study="Voolstra2023"
samples=$(cat ${prodir}/data/tara_samples.txt)
JOBFILE="$prodir/bash/jobs/fastqdump_${sample}.job"

for sample in $samples
do
echo "# /bin/sh" > $JOBFILE
#
echo "# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N fastqdump_${study}_${sample}
#$ -o $prodir/bash/jobs/fastqdump_${study}_${sample}.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/sratoolkit" >> $JOBFILE
echo "#" >> $JOBFILE
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE
#
echo "fastq-dump --split-files --outdir ${prodir}/data/raw/ --gzip $sample" >> $JOBFILE
#
echo "#" >> $JOBFILE
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done
