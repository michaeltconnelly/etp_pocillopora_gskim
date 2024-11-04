#!/bin/bash
#./bash/angsd_1DSFS_folded_reps_wrapper.sh
#purpose: calculate 1D SFS using ANGSD for a given list of populations
#input: ANGSD SAF files and list of populations pointing to files with lists of sample names

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)
echo "Begin 1D SFS calculation for all populations"

# script to create and submit separate jobs to calculate the 1D SFS for all species clusters/populations

for pop1 in $POPS; do
# create job file
echo "Creating job file for 1D SFS estimation of ${pop1}"
JOBFILE="${prodir}/bash/jobs/angsd_1DSFS_folded_reps_${pop1}.job"
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
#$ -N angsd_SFS_folded_${pop1}
#$ -o ${prodir}/bash/jobs/angsd_1DSFS_folded_reps_${pop1}.log
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
# NOTE: generate folded 1D-SFS with -fold 1 because ancestral state is not known
printf '/share/apps/bioinformatics/angsd/0.941/angsd/misc/realSFS ${angsddir}/safs/%s.reps.folded.saf.idx -fold 1 -P $NSLOTS > ${angsddir}/1dsfs/%s.reps.folded.ml \n' "$pop1" "$pop1" >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done

#for folded
#./misc/realSFS out.saf.idx -P 24 -fold 1 > out.sfs
#./misc/realSFS saf2theta out.saf.idx -outname out -sfs out.sfs -fold 1
#Estimate for every Chromosome/scaffold
#./misc/thetaStat do_stat out.thetas.idx
#Do a sliding window analysis based on the output from the make_bed command.
#./misc/thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz