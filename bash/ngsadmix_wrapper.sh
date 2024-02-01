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

# make output directory
mkdir -p ${prodir}/outputs/ngsadmix/${set}
for j in `seq 1 $K` ; do 
mkdir ${prodir}/outputs/ngsadmix/${set}/K${j}
done

# while loop to perform 5 iterations of NGSadmix analysis
i=1
while [ $i -le 5 ] ; do
# for loop to automate generation of scripts to perform NGSAdmix from 1:K 
for j in `seq 1 $K` ; do 
echo "NGSadmix for K=${j}, iteration #${i}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
echo "#$ -j y
#$ -N ngsadmix_${set}_K${j}_${i}.job
#$ -o ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
echo "module load bioinformatics/ngsadmix" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
#   input command for NGSadmix --> nested while loop runs for 10 iterations!
#
echo "NGSadmix -likes ${prodir}/outputs/angsd/${set}/${set}_ibs05.beagle.gz -K $j -P 16 -o ${prodir}/outputs/ngsadmix/${set}/K${j}/${set}_K${j}_${i} -minMaf 0.05" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
#
# submit job
qsub ${prodir}/bash/jobs/ngsadmix_${set}_K${j}_${i}.job
done
#
i=$(( $i + 1 ))
done

# NOTE input commands to organize log files in output for identifying the best K in R