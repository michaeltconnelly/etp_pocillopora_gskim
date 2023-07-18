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

# loop to automate generation of scripts to perform NGSAdmix from 2:K 
for j in `seq 2 $K`
do 
echo "NGSadmix for K=${j}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
echo "#$ -j y
#$ -N ngsadmix_${set}_K${j}.job
#$ -o ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
echo "module load bioinformatics/ngsadmix" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job

#   input command for NGSadmix
echo "NGSadmix -likes ${prodir}/outputs/angsd/${set}_ibs05.beagle.gz -K $j -P 16 -o ${prodir}/outputs/angsd/ngsadmix/${set} -minMaf 0.05" >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
# submit job
qsub ${prodir}/bash/jobs/ngsadmix_${set}_K${j}.job
done
