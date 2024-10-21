#!/bin/bash
#./bash/angsd_dxy_mafs_wrapper.sh
#purpose: calculate population maf posteriors using ANGSD for a given list of populations at specified sites
#input: list of populations pointing to files with lists of sample names, aligned bam files, site file 

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
# list of populations
POPFILE="$1"

# generate list of populations
POPS=$(cat ${prodir}/data/$POPFILE)

for set in $POPS; do 
# make bam file lists --> use only representative samples (top 5 for each species) to start
ls ${prodir}/outputs/alignments/pgra_himb/*md.rg.bam | grep -f ${prodir}/data/pops/${set} > ${angsddir}/bamfiles/${set}_all_bamfile.txt
# verify sample numbers are correct
echo "For ${set}, There are $(cat ${angsddir}/bamfiles/${set}_all_bamfile.txt | wc -l) samples in total"
# create job file
echo "Creating job file for population-specific MAF estimation of ${set}"
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
#$ -N angsd_dxy_maf_${set}
#$ -o ${prodir}/bash/jobs/angsd_dxy_maf_${set}.log
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
echo "#${set} SAF estimation" >> $JOBFILE
echo '# ----- SAF with ANGSD
# specify variable with bam file list for each population / species --> use only representative samples (top 5 for each species) to start
BAMS="${angsddir}/bamfiles/${set}_bamfile.txt" 

#FILTERS:
# restricting the MAF estimation to previously filtered sites (final_noclones_pgra_ibs05.sites.txt), using the -sites flag with a file corresponding to the recovered SNPs. This guarantees that sites with an allele fixed in one population are still included.
# -snp_pval 1e-5 : high confidence that the SNP is not just sequencing error 
# No -minInd, -minMaf filters because of previously filtered sites
# also adding filters against very badly non-HWE sites (such as, all calls are heterozygotes => lumped paralog situation) and sites with really bad strand bias:
FILTERS="-minMapQ 30 -minQ 30 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 "

#TODO: 
# -GL 1 : use samtools likelihood model
# -doMajorMinor 4 : uses the ref allele as major (-skipTriallelic 1 : skip tri-allelic sites)
# -doMaf 1 : estimate allele frequency w/ fixed major and minor alleles 
# -doCounts 1 : 
TODO="-doMajorMinor 4 -doMaf 1 -doCounts 1"

# set ancestral/reference genome fasta file - use reference for folded SFS spectra (ancestral state not known)
# NOTE: create folded SAF for 1-population analysis of neutrality test statistics
REF="/scratch/nmnh_corals/connellym/sequences/p_grandis_GCA_964027065.2/GCA_964027065.2_jaPocGran1.hap1.2_genomic.fna"' >> $JOBFILE

echo '# estimate the major allele frequency likelihood
angsd -sites ${angsddir}/final_noclones_pgra_ibs05.sites.txt -rf ${angsddir}/chrs_Pgra_HIMB.txt -b $BAMS -ref $REF -GL 1 -P $NSLOTS $TODO -out ${angsddir}/${set}

echo "DONE!"' >> $JOBFILE
# input job finished statment
echo '#
echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE $set
done