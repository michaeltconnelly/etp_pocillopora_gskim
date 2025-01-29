#!/bin/bash
#./bash/samtools_symbiont_C_D_fasta_jellyfish_wrapper.sh
#purpose: convert sorted, indexed symbiont_alignments to fasta read files with samtools, count kmers with jellyfish

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# Cladocopium
# making a list of sample names
set="C-dominated"
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the C-dominated samples to be processed:"
echo $samples

#loop to automate generation of scripts 
for sample in $samples
do \
JOBFILE="${prodir}/bash/jobs/${sample}_sym_C_samtools_split.job"
touch $JOBFILE
#
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_sym_C_samtools_jellyfish
#$ -o ${prodir}/bash/jobs/${sample}_sym_C_samtools_jellyfish.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/samtools" >> $JOBFILE
echo "module load tools/conda" >> $JOBFILE
echo "start-conda" >> $JOBFILE
echo "conda activate d2ssect" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

echo 'echo "Starting samtools fasta generation"' >> $JOBFILE

#   input command to produce fasta files for either Cladocopium or Durusdinium alignments
echo "samtools fasta ${prodir}/outputs/symbiont_alignments/${sample}_symC.sorted.bam \
-o ${prodir}/outputs/symbiont_kmers/C-dominated/${sample}_symC.fasta -@ 8"  >> $JOBFILE

#   input command for jellyfish kmer counting
echo "jellyfish count -m 21 -s 2000000000 ${prodir}/outputs/symbiont_kmers/C-dominated/${sample}_symC.fasta \
-o ${prodir}/outputs/symbiont_kmers/C-dominated/${sample}_symC.jf -t 8" >> $JOBFILE

echo 'echo '${sample}' successfully processed' >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
sleep 0.2
done

# Durusdinium
# making a list of sample names
set="D-dominated"
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the D-dominated samples to be processed:"
echo $samples

#loop to automate generation of scripts 
for sample in $samples
do \
JOBFILE="${prodir}/bash/jobs/${sample}_sym_D_samtools_split.job"
touch $JOBFILE
#
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > $JOBFILE
echo "# ----------------Parameters---------------------- #" >> $JOBFILE
echo "#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem" >> $JOBFILE
echo "#$ -j y
#$ -N ${sample}_sym_D_samtools_jellyfish
#$ -o ${prodir}/bash/jobs/${sample}_sym_D_samtools_jellyfish.log
#$ -m bea
#$ -M connellym@si.edu" >> $JOBFILE
#
echo "# ----------------Modules------------------------- #" >> $JOBFILE
echo "module load bioinformatics/samtools" >> $JOBFILE
echo "module load tools/conda" >> $JOBFILE
echo "start-conda" >> $JOBFILE
echo "conda activate d2ssect" >> $JOBFILE
#
echo "# ----------------Your Commands------------------- #" >> $JOBFILE
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> $JOBFILE
echo 'echo + NSLOTS = $NSLOTS' >> $JOBFILE

echo 'echo "Starting samtools fasta generation"' >> $JOBFILE

#   input command to produce fasta files for either Cladocopium or Durusdinium alignments
echo "samtools fasta ${prodir}/outputs/symbiont_alignments/${sample}_symD.sorted.bam \
-o ${prodir}/outputs/symbiont_kmers/D-dominated/${sample}_symD.fasta -@ 8"  >> $JOBFILE

#   input command for jellyfish kmer counting
echo "jellyfish count -m 21 -s 2000000000 ${prodir}/outputs/symbiont_kmers/D-dominated/${sample}_symD.fasta \
-o ${prodir}/outputs/symbiont_kmers/D-dominated/${sample}_symD.jf -t 8" >> $JOBFILE

echo 'echo '${sample}' successfully processed' >> $JOBFILE
#
echo 'echo = `date` job $JOB_NAME done' >> $JOBFILE
# submit job
qsub $JOBFILE
#
sleep 0.2
done