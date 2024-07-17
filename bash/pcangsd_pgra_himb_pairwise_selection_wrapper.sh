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
POPS=$(cat ${prodir}/data/$POPFILE)

# script to create and submit separate jobs for all pairwise comparisons
set -- $POPS
for pop1; do
shift
for pop2; do

# create job file
echo "Creating job file for PCangsd selection analysis between ${pop1} and ${pop2}"
JOBFILE="${prodir}/bash/jobs/pcangsd_pairwise_selection_${pop1}.${pop2}.job"
touch $JOBFILE

# create file with list of samples
SAMPLE_FILE="${prodir}/data/pairwise_comparisons/${pop1}_${pop2}_samples.txt"
cat ${prodir}/data/pops_ngsadmix/${pop1} ${prodir}/data/pops_ngsadmix/${pop2} > $SAMPLE_FILE

# input QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q mThC.q
#$ -l mres=8G,h_data=8G,h_vmem=8G
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
# assign job-specific variables
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"' >> $JOBFILE

# assign variables for pairwise comparison
printf 'set="%s_%s"\n' "$pop1" "$pop2" >> $JOBFILE

echo 'SAMPLE_FILE="${prodir}/data/pairwise_comparisons/${set}_samples.txt"' >> $JOBFILE

# input commands to create output directory and bam file list for each comparison
echo 'samples=$(cat $SAMPLE_FILE)

# creating an output directory for set output files
if [ ! -d "${angsddir}/selection/${set}" ]; then mkdir ${angsddir}/selection/${set}; fi

# assign variable for output directory
setdir="${angsddir}/selection/${set}"

# making a list of bam file paths
ls ${prodir}/outputs/alignments/pgra_himb/*md.rg.bam | grep -f $SAMPLE_FILE > ${setdir}/${set}_bamfile.txt

# verify samples are correct
echo "These are the samples to be processed: ${samples}, there are $(cat ${setdir}/${set}_bamfile.txt | wc -l) in total"' >> $JOBFILE

# input ANGSD commands
echo '# prepare ANGSD variables' >> $JOBFILE

# assign variable to sites and chrs files (sorted and indexed)
# include all sites, do not filter for LD
echo 'SITES="${angsddir}/AllSites_Pgra_HIMB.txt"
CHRS="${angsddir}/chrs_Pgra_HIMB.txt"' >> $JOBFILE

#echo 'SITES="${angsddir}/final_noclones/final_noclones_noLD_sorted_v3.sites.txt"
#CHRS="${angsddir}/final_noclones/final_noclones_noLD.SczhEnG_v3.txt"' >> $JOBFILE

# assign variable to bamfile
echo 'BAMS="${setdir}/${set}_bamfile.txt"' >> $JOBFILE

# assign variable to site filters 
echo 'FILTERS="-minMaf 0"' >> $JOBFILE

# assign variable to tasks
echo 'TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 8 -doPost 1 -doGlf 2"' >> $JOBFILE

# input ANGSD command to create BEAGLE file 
echo '# create the BEAGLE input file' >> $JOBFILE
echo 'angsd -sites $SITES -rf $CHRS -b $BAMS -GL 1 $FILTERS $TODO -P $NSLOTS -out ${setdir}/${set}_noLD_filtered' >> $JOBFILE

# input PCangsd commands
echo '# perform PCA-based selection analysis' >> $JOBFILE
echo 'pcangsd \
--beagle ${setdir}/${set}_noLD_filtered.beagle.gz \
-o ${setdir}/${set}_noLD \
--pcadapt \
--maf 0 \
--sites_save \
--threads $NSLOTS' >> $JOBFILE

# obtain site coordinates from finished maf.gz file
echo 'zcat ${setdir}/${set}_noLD_filtered.mafs.gz | cut -f 1,2 | tail -n +2 | sed 's/OZ//g' | sort -t$'\t' -k 1,1n -k 2,2n | sed 's/^/OZ/g' > ${setdir}/${set}_noLD.sites.txt'  >> $JOBFILE

# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE

# submit job
qsub $JOBFILE
done
done