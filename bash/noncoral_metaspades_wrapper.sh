#!/bin/bash
#./bash/noncoral_metaspades_wrapper.sh
#purpose: spades assembly

# specify variable containing sequence file prefixes and directory paths to shorten paths 
# scratch directory 
mcs="/scratch/nmnh_corals/connellym" 
# project root working directory
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim" 

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt) # this can be any file that has a list of sample file prefixes

#lets me know which files are being processed
echo "These are the samples to be assembled:"
echo $samples

# loop to automate generation of scripts to direct sequence file trimming
# edit output job script location within your project root directory
JOBFILE="${prodir}/bash/jobs/${sample}_noncoral_metaspades.job"
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=384G,h_data=24G,h_vmem=24G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_noncoral_metaspades
#$ -o ${prodir}/bash/jobs/${sample}_noncoral_metaspades.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/spades" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

#   input command for metaspades assembly of non-coral reads
echo "spades.py --meta \
-o ${prodir}/outputs/noncoral_spades/${sample} \
--pe1-1 ${prodir}/data/unmapped/${sample}_Unmapped_R1_PE_trimmed.fastq.gz \
--pe1-2 ${prodir}/data/unmapped/${sample}_Unmapped_R2_PE_trimmed.fastq.gz"  >> $JOBFILE
#
echo 'echo '${sample}' successfully assembled' >> "${prodir}"/bash/jobs/${sample}_spades.job
#
#   move contigs, cleanup intermediate files
echo "cp ${prodir}/outputs/noncoral_spades/${sample}/contigs.fasta ${prodir}/outputs/noncoral_spades/contigs/${sample}_contigs.fasta" >> $JOBFILE

echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
done