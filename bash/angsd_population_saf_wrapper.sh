#!/bin/bash
#./bash/angsd_population_saf_wrapper.sh
#purpose: calculate population site allele frequency posteriors using ANGSD for a given list of populations
#input: list of populations pointing to files with lists of sample names, aligned bam files

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/pops/$POPFILE)

# make bam file lists
for set in $POPS; do 
ls ${prodir}/outputs/alignments/*md.rg.bam | grep -f ${prodir}/data/pops/${set} > ${angsddir}/bamfiles/${set}_bamfile.txt
# verify sample numbers are correct
echo "For ${set}, There are $(cat ${angsddir}/bamfiles/${set}_bamfile.txt | wc -l) samples in total"

# create job file
echo "Creating job file for SAF estimation of ${set}"
JOBFILE="${prodir}/bash/jobs/angsd_saf_${set}.job"
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
#$ -N angsd_saf_${set}
#$ -o ${prodir}/bash/jobs/angsd_saf_${set}.log
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
set="$1"
#' >> $JOBFILE

# input ANGSD commands
echo "${set} SAF estimation" >> $JOBFILE
echo '# ----- SAF with ANGSD
# specify variable with bam file list for each population / species
BAMS="${angsddir}/bamfiles/${set}_bamfile.txt"
# filters are not needed since we are restricting the SAF estimation to previously filtered sites (AllSites)
# -doSaf 1: perform multisample GL estimation
TODO=" -doSaf 1"
# set ancestral reference genome fasta file
#ANC="/scratch/nmnh_corals/connellym/sequences/pdam/pdam_genome.fasta"
ANC="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/seqs/peffusa_reference_ancestral.fa"' >> $JOBFILE

echo '# estimate the site allele frequency likelihood
angsd -sites ${angsddir}/AllSites.txt -rf ${angsddir}/SczhEnG.txt -b $BAMS -anc $ANC -GL 1 -P $NSLOTS $TODO -out ${angsddir}/${set} 
echo "DONE!"' >> $JOBFILE
#qsub $JOBFILE $set
done