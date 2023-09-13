#!/bin/bash

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
study="$1"

echo "# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=16G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N symalignqc_${study}
#$ -o ${prodir}/bash/jobs/symalignqc_${study}.log
#$ -m bea
#$ -M connellym@si.edu" > $prodir/bash/jobs/symalignqc_${study}.job
#
echo "#" >> $prodir/bash/jobs/symalignqc_${study}.job
echo "# ----------------Modules------------------------- #" >> $prodir/bash/jobs/symalignqc_${study}.job
echo "module load bioinformatics/samtools" >> $prodir/bash/jobs/symalignqc_${study}.job
echo "#
# ----------------Your Commands------------------- #
#" >> $prodir/bash/jobs/symalignqc_${study}.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $prodir/bash/jobs/symalignqc_${study}.job
echo 'echo + NSLOTS = $NSLOTS' >> $prodir/bash/jobs/symalignqc_${study}.job
echo "#" >> $prodir/bash/jobs/symalignqc_${study}.job
echo 'prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"' >> $prodir/bash/jobs/symalignqc_${study}.job
# produce samtools flagstat reports for alignment QC
echo 'for i in ${prodir}/outputs/symbiont_alignments/*.sam; do name=$(basename $i | cut -d . -f 1); echo $name; samtools flagstat $i > ${prodir}/outputs/QCs/symbiont_flagstats/${name}_symbiont_flagstat.txt; done' >> $prodir/bash/jobs/symalignqc_${study}.job

echo "echo = `date` job $JOB_NAME done" >> $prodir/bash/jobs/symsymalignqc_${study}.job

qsub $prodir/bash/jobs/symalignqc_${study}.job

