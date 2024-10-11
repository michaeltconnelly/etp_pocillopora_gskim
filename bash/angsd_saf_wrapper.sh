#!/bin/bash
#./bash/angsd_saf_wrapper.sh
#purpose: calculate population site allele frequency posteriors using ANGSD for a given list of populations
#input: list of populations pointing to files with lists of sample names, aligned bam files

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)

for set in $POPS; do 
# make bam file lists
ls ${prodir}/outputs/alignments/pgra_himb/*md.rg.bam | grep -f ${prodir}/data/pops_ngsadmix/${set} > ${angsddir}/bamfiles/${set}_bamfile.txt
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
module load bioinformatics/angsd
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
# set ancestral/reference genome fasta file - use reference for folded SFS spectra (ancestral state not known)
ANC="/scratch/nmnh_corals/connellym/sequences/p_grandis_GCA_964027065.2/GCA_964027065.2_jaPocGran1.hap1.2_genomic.fna"' >> $JOBFILE
echo '# estimate the site allele frequency likelihood
angsd -sites ${angsddir}/AllSites_Pgra_HIMB.txt -rf ${angsddir}/chrs_Pgra_HIMB.txt -b $BAMS -anc $ANC -GL 1 -P $NSLOTS $TODO -out ${angsddir}/${set} 
echo "DONE!"' >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE $set
done