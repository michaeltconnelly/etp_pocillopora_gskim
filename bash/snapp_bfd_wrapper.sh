#!/bin/bash
#./bash/snapp_bfd_wrapper.sh
#purpose: submit jobs for SNAPP runs with different species hypotheses for BFD*
#input: properly configured XML files 

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
bfddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/snapp/BFD"
# specify which run xml file to use (runA, runB, etc.)
run=$1

# create job file
echo "Creating job for SNAPP ${run}"
JOBFILE="${bfddir}/${run}/snapp_${run}.job"
touch $JOBFILE
# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 48
#$ -q lThC.q
#$ -l mres=96G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N SNAPP_path_${run}
#$ -o ${bfddir}/${run}/snapp_path_${run}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/beast/2.7.6-beagle
# ----------------Your Commands------------------- #
#" > $JOBFILE

# input job-specific variables
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
bfddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/snapp/BFD"
run="$1"
#' >> $JOBFILE

# input SNAPP commands
# note - call the qsub command from the intended output directory
echo 'cd ${bfddir}/${run}' >> $JOBFILE
echo '# start SNAPP job
beast -threads $NSLOTS ${run}.xml
echo "DONE!"' >> $JOBFILE

# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE

# submit job
qsub $JOBFILE $run