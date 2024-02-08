#!/bin/bash
#./bash/angsd_pairwise_fst_wrapper.sh
#purpose: calculate pairwise fst values using ANGSD for a given list of populations
#input: ANGSD SAF and pairwise SFS files and list of populations pointing to files with lists of sample names

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/pops/$POPFILE)
# calculate genome-wide per-site pairwise Fst
# script to create and submit separate jobs for all pairwise comparisons

set -- $POPS
for pop1; do
shift
for pop2; do
# create job file
echo "Creating job file for pairwise Fst comparison of ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/angsd_pairwise_fst_${pop1}.${pop2}.job"
touch $JOBFILE
# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -q lThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N angsd_pairwise_fst_${pop1}.${pop2}
#$ -o ${prodir}/bash/jobs/angsd_pairwise_fst_${pop1}.${pop2}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bio/angsd/0.940
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
printf 'realSFS fst index ${angsddir}/%s.saf.idx ${angsddir}/%s.saf.idx -sfs ${angsddir}/%s.%s.ml -fstout ${angsddir}/%s_%s \n' "$pop1" "$pop2" "$pop1" "$pop2" "$pop1" "$pop2" >> $JOBFILE
echo '# get the global Fst estimate' >> $JOBFILE
printf 'realSFS fst stats ${angsddir}/%s_%s.fst.idx > ${angsddir}/%s_%s_fst_results.txt \n' "$pop1" "$pop2" "$pop1" "$pop2" >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done
done