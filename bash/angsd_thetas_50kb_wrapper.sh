#!/bin/bash
#./bash/angsd_thetas_50kb_wrapper.sh
#purpose: calculate window-based thetas from folded 1D-SFS using ANGSD for a given list of populations
#input: ANGSD SAF files and list of populations pointing to files with lists of sample names

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)
echo "Begin calculation of thetas from 1D-SFS for all populations"

# script to create and submit separate jobs to calculate the 2D SFS for all pairwise comparisons

for pop1 in $POPS; do
# create job file
echo "Creating job file for theta estimation of ${pop1}"
JOBFILE="${prodir}/bash/jobs/angsd_thetas_${pop1}.job"
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
#$ -N angsd_thetas_50kb_${pop1}
#$ -o ${prodir}/bash/jobs/angsd_thetas_50kb_${pop1}.log
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
# NOTE: calculate per-site thetas using folded 1D-SFS because ancestral state is not known
# estimate theta in 50kb windows
printf '/share/apps/bioinformatics/angsd/0.941/angsd/misc/thetaStat do_stat ${angsddir}/thetas/%s.thetas.idx -win 50000 -step 50000 -type 2 -outnames ${angsddir}/thetas/windows_50kb/%s.thetas.50kb.non-overlap \n' "$pop1" "$pop1" >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done