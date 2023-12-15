#!/bin/bash

# generate samtools flagstat reports for comparing read alignment rates to different Pocillopora reference genomes
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
set="$1"
genome="$2"

echo "# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=16G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N alignqc_${set}_${genome}
#$ -o ${prodir}/bash/jobs/alignqc_${set}_${genome}.log
#$ -m bea
#$ -M connellym@si.edu" > $prodir/bash/jobs/alignqc_${set}_${genome}.job
#
echo "#" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo "# ----------------Modules------------------------- #" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo "module load bioinformatics/samtools" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo "#
# ----------------Your Commands------------------- #
#" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo 'echo + NSLOTS = $NSLOTS' >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo "#" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
echo 'prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"' >> $prodir/bash/jobs/alignqc_${set}_${genome}.job
# produce samtools flagstat reports for alignment QC
echo 'for i in ${prodir}/outputs/alignments/${genome}/*.sam; do name=$(basename $i | cut -d . -f 1); echo $name; samtools flagstat $i > ${prodir}/outputs/QCs/flagstats/${genome}/${name}_${genome}_flagstat.txt; done' >> $prodir/bash/jobs/alignqc_${set}_${genome}.job

echo "echo = `date` job $JOB_NAME done" >> $prodir/bash/jobs/alignqc_${set}_${genome}.job

qsub $prodir/bash/jobs/alignqc_${set}_${genome}.job

