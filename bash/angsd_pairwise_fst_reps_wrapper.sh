#!/bin/bash
#./bash/angsd_pairwise_fst_reps_wrapper.sh
#purpose: calculate pairwise fst values using ANGSD for a given list of populations
#input: ANGSD SAF and pairwise SFS files and list of populations pointing to files with lists of sample names

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)
# calculate genome-wide per-site pairwise Fst
# script to create and submit separate jobs for all pairwise comparisons

set -- $POPS
for pop1; do
shift
for pop2; do
# create job file
echo "Creating job file for pairwise Fst comparison of ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/angsd_pairwise_fst_reps_${pop1}.${pop2}.job"
touch $JOBFILE
# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -q lThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N angsd_pairwise_fst_reps_${pop1}.${pop2}
#$ -o ${prodir}/bash/jobs/angsd_pairwise_fst_reps_${pop1}.${pop2}.log
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
echo '# prepare the Fst' >> $JOBFILE
printf '/share/apps/bioinformatics/angsd/0.941/angsd/misc/realSFS fst index ${angsddir}/safs/%s.reps.folded.saf.idx ${angsddir}/safs/%s.reps.folded.saf.idx -sfs ${angsddir}/2dsfs/%s.%s.reps.folded.ml -fold 1 -fstout ${angsddir}/fst/%s_%s.reps \n' "$pop1" "$pop2" "$pop1" "$pop2" "$pop1" "$pop2" >> $JOBFILE
echo '# get the global Fst estimate' >> $JOBFILE
printf '/share/apps/bioinformatics/angsd/0.941/angsd/misc/realSFS fst stats ${angsddir}/fst/%s_%s.reps.fst.idx > ${angsddir}/fst/%s_%s_reps_fst_results.txt \n' "$pop1" "$pop2" "$pop1" "$pop2" >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
done
done