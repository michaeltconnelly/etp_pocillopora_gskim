#!/bin/bash

# symbiont barcodes to try: mtcob, mtcox1, ITS1, ITS2, LSU, cp23S, psbAncr
# Cladocopium goreaui, C. latusorum, C. pacificum, Durusdinium glynni, D. trenchii


for infile in ${prodir}/outputs/noncoral_spades/contigs/*_contigs.fasta
do
     # Make BLAST database from non-coral assemblies
     base=$(basename ${infile} _contigs.fasta)
     makeblastdb -in ${infile} -out ${prodir}/outputs/noncoral_spades/blastdbs/${base} -dbtype nucl -parse_seqids 
     # mt COB genes
     # C. latusorum
     blastn -query ${prodir}/data/seqs/C_lat_mtcob_MW713374.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_C_lat_mtcob
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_C_lat_mtcob | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_C_lat_mtcob.sorted
     # D. glynni
     blastn -query ${prodir}/data/seqs/D_gly_mtcob_KY131780.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_D_gly_mtcob
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_D_gly_mtcob | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/mtcob_blast_hits/${base}_D_gly_mtcob.sorted
done


# Put top hits from all samples into one text file
# python blast_sum.py