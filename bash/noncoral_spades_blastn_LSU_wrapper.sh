#!/bin/bash

# symbiont barcodes to try: LSU, mtcox1, ITS1, ITS2, LSU, cp23S, psbAncr
# Cladocopium goreaui, C. latusorum, C. pacificum, Durusdinium glynni, D. trenchii

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

for infile in $(ls ${prodir}/outputs/noncoral_spades/contigs/*_contigs.fasta)
do
     # Use previously made BLAST databases from non-coral assemblies 
     base=$(basename ${infile} _contigs.fasta)
     # LSU genes
     # C. latusorum
     blastn -query ${prodir}/data/seqs/C_lat_LSU_MW711730.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_lat_LSU
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_lat_LSU | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_lat_LSU.sorted
     # C. pacificum 
     blastn -query ${prodir}/data/seqs/C_pac_LSU_MW711733.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_pac_LSU
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_pac_LSU | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_C_pac_LSU.sorted
     # D. glynni
     blastn -query ${prodir}/data/seqs/D_gly_LSU_KY131781.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_D_gly_LSU
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_D_gly_LSU | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/LSU_blast_hits/${base}_D_gly_LSU.sorted
done


# Put top hits from all samples into one text file
# python blast_sum.py