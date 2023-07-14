#!/bin/bash
# purpose: after running SPAdes, extract PocHistone3 sequence from assembled contigs using BLAST

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
samples=$(cat ${prodir}/data/${set}_samples.txt)

# 
for sample in $samples
do
# make BLAST database to search for histone3 sequences
makeblastdb -in ${prodir}/outputs/spades/${sample}*/contigs.fasta -dbtype 'nucl' -title ${prodir}/outputs/spades/${sample}*/contigs -parse_seqids

# perform BLAST search
blastn -query ${prodir}/data/seqs/Pmeandrina_Hist3_MG587097.fasta -db ${prodir}/outputs/spades/${sample}*/contigs.fasta -out ${prodir}/outputs/pochist/blast_results/${sample}_PocHist3.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -max_target_seqs 1

# extract best BLAST hit subject sequence as fasta
cat ${prodir}/outputs/pochist/blast_results/${sample}_PocHist3.tsv | cut -d$'\t' -f 2,13 | sed "s/^/>${sample}_/g" | sed 's/\t/\n/g' > ${prodir}/outputs/pochist/${sample}_hist3.fasta

# collate results
cat ${prodir}/outputs/pochist/${sample}_hist3.fasta >> ${prodir}/outputs/pochist/PocHist3_All.fasta
done