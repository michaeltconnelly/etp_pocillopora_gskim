#!/bin/bash

# specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# setup loop to create sharkmer jobs
set=$1
cat ${prodir}/data/${set}_samples.txt | while read SAMPLE 

do 
# create job file 
JOBFILE=${prodir}/bash/jobs/${SAMPLE}_sharkmer.job
# add QSUB commands
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 2
#$ -q sThM.q
#$ -l mres=80G,h_data=40G,h_vmem=40G,himem
#$ -cwd
#$ -j y
#$ -N sharkmer_${SAMPLE}
#$ -o ${prodir}/outputs/sharkmer/${SAMPLE}/sharkmer_${SAMPLE}.log
#" > $JOBFILE
# add modules 
echo '# ----------------Modules------------------------- #
  module load bio/sharkmer/0.2.0_93ee045
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE

echo "gunzip ${prodir}/data/trimmed/${SAMPLE}_R1_PE_trimmed.fastq.gz
gunzip ${prodir}/data/trimmed/${SAMPLE}_R2_PE_trimmed.fastq.gz" >> $JOBFILE

echo "sharkmer \
 --max-reads 6000000 -k 33 \
 -s ${SAMPLE} -o ${prodir}/outputs/sharkmer/${SAMPLE}/ \
  --pcr "TTTGGGSATTCGTTTAGCAG,SCCAATATGTTAAACASCATGTCA,1500,mtORF,coverage=20,mismatches=1,trim=24" \
  --pcr "ATTCAGTCTCACTCACTCACTCAC,TATCTTCGAACAGACCCACCAAAT,1000,PocHistone,coverage=20,mismatches=1,trim=24" \
 ${prodir}/data/trimmed/${SAMPLE}_R1_PE_trimmed.fastq ${prodir}/data/trimmed/${SAMPLE}_R2_PE_trimmed.fastq" >> $JOBFILE

echo "gzip ${prodir}/data/trimmed/${SAMPLE}_R1_PE_trimmed.fastq
gzip ${prodir}/data/trimmed/${SAMPLE}_R2_PE_trimmed.fastq" >> $JOBFILE

echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
done 
#
#TO RUN JOB - 
# All reads files should be unzipped (use gunzip command to do this)
# Include unique sample IDs in file names of read files = {SAMPLE}. 
# Create a text file called "sharkmer_samples.txt" with each unique sample name in a
# column (i.e., one sample name per line).
# Modify the input read files in script (e.g., ./${SAMPLE}_R1_PE_trimmed.fastq) as needed
# to match your file names for forward and reverse reads (usually contain R1 and R2).
# Run job file in same directory as read files and sharkmer_samples.txt file.

# NOTES -
# To add or remove primer sequences, edit the --pcr commands in script. Include forward
# and reverse primer sequences, estimated length (overestimate this length), and primer
# pair name.

