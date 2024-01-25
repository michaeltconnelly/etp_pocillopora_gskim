#!/bin/bash
#./bash/vcftools_pairwise_fst_wrapper.sh
#purpose: calculate pairwise fst values using vcftools for a given list of populations
#input: strict-filtered, hard-called VCF file and list of populations pointing to files with lists of sample names

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
# list of populations
VCF_NAME="$1"
POPFILE="$2"

touch ${prodir}/bash/jobs/vcftools_pairwise_fst_${POPFILE}.job

echo "Preparing script for ${POPS}"
#   input QSUB commands
echo '#!/bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -q lThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N vcftools_pairwise_fst
#$ -o vcftools_pairwise_fst.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/vcftools
#
# ----------------Your Commands------------------- #
#
#' > ${prodir}/bash/jobs/vcftools_pairwise_fst_${POPFILE}.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
VCF_NAME="$1"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
VCF_IN="${prodir}/outputs/${VCF_NAME}.vcf.gz"
#
#' >> ${prodir}/bash/jobs/vcftools_pairwise_fst_${POPFILE}.job
#
# generate list of populations
POPS=$(cat ${prodir}/data/pops/$POPFILE)
# calculate genome-wide per-site pairwise Fst
# script to create commands for all pairwise comparisons
set -- $POPS
for a; do
    shift
    for b; do
        printf 'vcftools --gzvcf $VCF_IN --weir-fst-pop ${prodir}/data/pops/%s --weir-fst-pop ${prodir}/data/pops/%s --out ${prodir}/outputs/%s_%s\n' "$a" "$b" "$a" "$b" >> ${prodir}/bash/jobs/vcftools_pairwise_fst_${POPFILE}.job
    done
done
#

echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/vcftools_pairwise_fst_${POPFILE}.job