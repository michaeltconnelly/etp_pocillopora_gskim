#!/bin/bash
#./bash/angsd_pairwise_dxy_wrapper.sh
#purpose: calculate pairwise dxy using ANGSD for a given list of populations
#input: ANGSD population-specific MAF files and list of populations pointing to files with lists of sample names

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)
echo "Begin dxy calculation for all combinations"

# script to create and submit separate jobs to calculate the 2D SFS for all pairwise comparisons

set -- $POPS
for pop1; do
shift
for pop2; do
# create job file
echo "Creating job file for pairwise dxy estimation of ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/angsd_dxy_${pop1}.${pop2}.job"
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
#$ -o ${prodir}/bash/jobs/angsd_dxy_${pop1}.${pop2}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/angsd
module load tools/R
#
# ----------------Your Commands------------------- #
#" > $JOBFILE
# input job-specific variables
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
#' >> $JOBFILE
# run analysis from ANGSD output directory
echo 'cd $angsddir' >> $JOBFILE
# input ANGSD and auxiliary script commands
# script adapted from instructions 
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/calcDxy.R
# https://github.com/ksamuk/pixy_analysis/blob/main/data_generation/angsd/ANGSD_output/README_ANGSDforPixyPaper.txt
# STEPS:
# 1) Run ANGSD with all populations with -SNP_pval and -skipTriallelic flags. Use the ref allele as the major allele ("-doMajorMinor 4"). Discard genotypes with Q score below 30 or depth below 10.
# Completed in angsd_ibs.job, 5,587,480 study-wide SNPs after filtering
# 2) Re-run ANGSD per population, using the -sites flag with a file corresponding to the recovered SNPs. This guarantees that sites with an allele fixed in one population are still included. 
# Completed in angsd_dxy_mafs_wrapper.sh
# 3) Gunzip the resulting mafs files and run calcDxy.R script
printf 'Rscript ${prodir}/R/calcDxy.R -p %s.mafs -q %s.mafs \n' "$pop1" "$pop2" >> $JOBFILE
# This step makes “Dxy_persite.txt" and outputs a global Dxy
echo "mv ${angsddir}/Dxy_persite.txt ${angsddir}/dxy/${pop1}_${pop2}_Dxy_persite.txt" >> $JOBFILE

# 4) Run custom perl script to average the persite output over 50kb windows.
echo "perl ${prodir}/bash/angsd_dxy_step3_processoutputoverwindows_04242020.pl ${angsddir}/dxy/${pop1}_${pop2}_Dxy_persite.txt 50000" >> $JOBFILE
# The output is chrX_angsd_DxySummary.txt, which includes an avg Dxy per window using the total # of sites as the denominator, and an estimate using the # of sites with data as the denominator.
echo "mv ${angsddir}/chrX_angsd_DxySummary.txt ${angsddir}/dxy/${pop1}_${pop2}_chrX_angsd_DxySummary.txt" >> $JOBFILE

# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
qsub $JOBFILE
done
done