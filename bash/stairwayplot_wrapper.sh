#!/bin/bash
#./bash/stairwayplot_wrapper.sh
#purpose: prepare a job file to run stairwayplot on Hydra
#input: ANGSD folded 1D SFS files and stairwayplot v2 blueprint file with configuration options

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
stairwayplot="/scratch/nmnh_corals/connellym/programs/stairway_plot_v2.1.1/stairway_plot_es"
pop1="$1"

java -cp $stairwayplot Stairbuilder ${prodir}/bash/stairwayplot_blueprints/${pop1}_fold.blueprint

echo "Creating job file for stairway plot using 1D SFS of ${pop1}"
JOBFILE="${prodir}/bash/stairwayplot_blueprints/stairwayplot_folded_${pop1}.job"
touch $JOBFILE

# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q lThM.q
#$ -l mres=12G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N stairwayplot_folded_${pop1}
#$ -o ${prodir}/bash/jobs/stairwayplot_folded_${pop1}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
#
# ----------------Your Commands------------------- #
#" > $JOBFILE

# input job-specific variables
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
#' >> $JOBFILE

# input stairwayplot commands
cat ${prodir}/bash/stairwayplot_blueprints/${pop1}_fold.blueprint.sh >> $JOBFILE

# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE

# submit job file
# qsub $JOBFILE