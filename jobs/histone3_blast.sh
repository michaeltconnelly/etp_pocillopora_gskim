#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N histone3_blastloop
#$ -o histone3_blastloop.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/blast
#
# ----------------Your Commands------------------- #
# purpose: after running SPAdes, extract PocHistone3 sequence from assembled contigs using BLAST
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
samples=$(cat ${prodir}/data/${set}_samples.txt)

# make output directory structure
mkdir -p ${prodir}/outputs/pochist/blast_results/ ${prodir}/outputs/pochist/fastas/

# 
for sample in $samples
do
# assign variable for sample spades assembly directory
spades_dir=${prodir}/outputs/spades/${sample}*

# add sample ID string to contig fasta headers
sed "s/^>/>${sample}_/g" ${spades_dir}/contigs.fasta | sed 's/ODE_//g' > ${prodir}/outputs/spades/${sample}_contigs.fasta
done

for sample in $samples
do
# make BLAST database to search for histone3 sequences
makeblastdb -in ${prodir}/outputs/spades/contigs/${sample}_contigs.fasta -out ${prodir}/outputs/spades/blastdbs/${sample} -dbtype 'nucl' -title ${sample} -parse_seqids

# perform BLAST search
blastn -query ${prodir}/data/seqs/Pmeandrina_Hist3_MG587097.fasta -db ${prodir}/outputs/spades/blastdbs/${sample} -out ${prodir}/outputs/pochist/blast_results/${sample}_PocHist3.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -max_target_seqs 1

# extract best BLAST hit subject sequence as fasta
cat ${prodir}/outputs/pochist/blast_results/${sample}_PocHist3.tsv | cut -d$'\t' -f 2,13 | sed "s/^/>${sample}_/g" | sed 's/\t/\n/g' > ${prodir}/outputs/pochist/fastas/${sample}_hist3.fasta

# collate results
cat ${prodir}/outputs/pochist/fastas/${sample}_hist3.fasta >> ${prodir}/outputs/pochist/PocHist3_All.fasta
done