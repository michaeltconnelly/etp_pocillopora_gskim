#!/bin/bash

# specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# setup loop to create kraken jobs
set=$1
DBNAME="${prodir}/data/krakendb/k2_standard_16gb_20240605"
#DBNAME="${prodir}/data/krakendb/k2_pluspf_20240605"

cat ${prodir}/data/${set}_samples.txt | while read SAMPLE 
do 
# create job file 
JOBFILE=${prodir}/bash/jobs/${SAMPLE}_kraken_spades.job

# make sample output directory
OUTDIR="${prodir}/outputs/kraken/${SAMPLE}"

mkdir $OUTDIR

# add QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 12
#$ -q sThM.q
#$ -l mres=144G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N kraken_classify_${SAMPLE}_spades
#$ -o ${prodir}/outputs/kraken/${SAMPLE}/kraken_${SAMPLE}_spades.log
#$ -m bea
#$ -M connellym@si.edu
#" > $JOBFILE

# add modules 
echo '# ----------------Modules------------------------- #
module load bioinformatics/kraken
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE

echo "kraken2 --threads 12 --db $DBNAME \
--output ${OUTDIR}/${SAMPLE}.tsv \
--classified-out ${OUTDIR}/${SAMPLE}_classified.fa \
--unclassified-out ${OUTDIR}/${SAMPLE}_unclassified.fa \
--report ${OUTDIR}/${SAMPLE}_report \
--use-names \
${prodir}/outputs/spades/contigs/${set}/${SAMPLE}_contigs.fasta" >> $JOBFILE


echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE

done 