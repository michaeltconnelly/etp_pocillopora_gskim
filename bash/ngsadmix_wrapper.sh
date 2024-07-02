#!/bin/bash
#./bash/ngsadmix_wrapper.sh
# purpose: wrapper script to submit jobs that perform NGSadmix on multiple values of K ancestral populations

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
# set a value for max number of K populations
K=$2

# make output directory structure
mkdir -p ${prodir}/outputs/ngsadmix/${set}
for j in `seq 1 $K` ; do 
mkdir ${prodir}/outputs/ngsadmix/${set}/K${j}
done

# while loop to perform 10 iterations of NGSadmix analysis
i=1
while [ $i -le 10 ] ; do
# for loop to automate generation of scripts to perform NGSadmix from 1:K 
for j in `seq 1 $K` ; do 
echo "NGSadmix for K=${j}, iteration #${i}"
# create job file
JOBFILE="${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job"
touch $JOBFILE
# input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem
#$ -j y
#$ -N ngsadmix_${set}_K${j}_${i}.job
#$ -o ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.log
#$ -m bea
#$ -M connellym@si.edu
# ----------------Modules------------------------- #
module load bio/ngsadmix/32
# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS' >> $JOBFILE
#
# input command for NGSadmix 
echo "NGSadmix -likes ${prodir}/outputs/angsd/${set}/${set}_noLD_pgra_himb_filtered.beagle.gz \
-K $j -P 16 \
-o ${prodir}/outputs/ngsadmix/${set}/K${j}/${set}_noLD_pgra_himb_K${j}_${i} \
-minMaf 0.05" >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
#
# submit job
qsub $JOBFILE
done
#
i=$(( $i + 1 ))
done

# NOTE input commands to organize log files in output for identifying the best K in R