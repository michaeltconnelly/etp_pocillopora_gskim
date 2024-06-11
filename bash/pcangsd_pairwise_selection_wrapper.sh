#!/bin/bash
#./bash/pcangsd_pairwise_selection_wrapper.sh
#purpose: perform PCA-based selection analyses between pairs of populations
#input: lists of sample names

# specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"

# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/pops_ngsadmix/$POPFILE)

# script to create and submit separate jobs for all pairwise comparisons
set -- $POPS
for pop1; do
shift
for pop2; do
# create job file
echo "Creating job file for PCangsd selection analysis between ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/pcangsd_pairwise_selection_${pop1}.${pop2}.job"
touch $JOBFILE

# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q lThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N pcangsd_selection_${pop1}.${pop2}
#$ -o ${prodir}/bash/jobs/pcangsd_selection_${pop1}.${pop2}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/angsd
module load bioinformatics/pcangsd
#
# ----------------Your Commands------------------- #
#" > $JOBFILE

# input job-specific variables
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
#' >> $JOBFILE

# assign variable for pairwise comparison
printf 'set="%s_%s"\n' "$pop1" "$pop2" >> $JOBFILE

# input commands to create output directory and bam file list for each comparison
echo 'samples=$(cat ${prodir}/data/${set}_samples.txt)

# creating an output directory for set output files
if [ ! -d "${angsddir}/${set}" ]; then mkdir ${angsddir}/${set}; fi

# assign variable for output directory
setdir="${angsddir}/${set}"

# making a list of bam file paths
ls ${prodir}/outputs/alignments/*md.rg.bam | grep -f ${prodir}/data/${set}_samples.txt > ${setdir}/${set}_bamfile.txt

# verify samples are correct
echo " These are the samples to be processed: ${samples}, there are $(cat ${setdir}/${set}_bamfile.txt | wc -l) in total"' >> $JOBFILE

# input ANGSD commands
echo '# prepare ANGSD variables' >> $JOBFILE

# assign variable to sites and chrs files (sorted and indexed)
echo 'SITES=${angsddir}/final_noclones/final_noclones_noLD_sorted_v3.sites.txt
CHRS=${angsddir}/final_noclones/final_noclones_noLD.SczhEnG_v3.txt"' >> $JOBFILE

# assign variable to bamfile
echo 'BAMS=${setdir}/${set}_bamfile.txt' >> $JOBFILE

# assign variable to site filters 
echo 'FILTERS="-snp_pval 1e-5 -minMaf 0.05"' >> $JOBFILE

# assign variable to tasks
echo 'TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 8 -doPost 1 -doGlf 2"' >> $JOBFILE

# input ANGSD command to create BEAGLE file 
echo '# create the BEAGLE input file' >> $JOBFILE
echo 'angsd -sites $SITES -rf $CHRS -b $BAMS -GL 1 $FILTERS $TODO -P $NSLOTS -out ${setdir}/${set}_noLD_filtered' >> $JOBFILE

# input PCangsd commands
echo '# perform PCA-based selection analysis' >> $JOBFILE
echo 'pcangsd \
-beagle ${angsddir}/${set}/${set}_noLD_filtered.beagle.gz \
-o ${prodir}/outputs/angsd/selection/${set}_noLD \
--pcadapt \
--sites_save \
--threads $NSLOTS' >> $JOBFILE

# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE

# submit job
# qsub $JOBFILE
done
done