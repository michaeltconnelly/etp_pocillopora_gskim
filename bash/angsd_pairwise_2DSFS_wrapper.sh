#!/bin/bash
#./bash/angsd_pairwise_2DSFS_wrapper.sh
#purpose: calculate pairwise 2D SFS using ANGSD for a given list of populations
#input: ANGSD SAF files and list of populations pointing to files with lists of sample names

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/pops_ngsadmix/$POPFILE)
echo "Begin 2D SFS calculation for all combinations"

# script to create and submit separate jobs to calculate the 2D SFS for all pairwise comparisons

set -- $POPS
for pop1; do
shift
for pop2; do
# create job file
echo "Creating job file for pairwise 2D SFS estimation of ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/angsd_2DSFS_${pop1}.${pop2}.job"
touch $JOBFILE
# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N angsd_2DSFS_${pop1}.${pop2}
#$ -o ${prodir}/bash/jobs/angsd_2DSFS_${pop1}.${pop2}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/angsd
#
# ----------------Your Commands------------------- #
#" > $JOBFILE
# input job-specific variables
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
#' >> $JOBFILE
# input ANGSD commands
printf 'realSFS ${angsddir}/%s.saf.idx ${angsddir}/%s.saf.idx -P $NSLOTS > ${angsddir}/%s.%s.ml \n' "$pop1" "$pop2" "$pop1" "$pop2" >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done
done