#!/bin/bash

# symbiont barcode: ITS2
# C. latusorum, C. pacificum, D.glynni

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

for infile in $(ls ${prodir}/outputs/noncoral_spades/contigs/*_contigs.fasta)
do
     # Use previously made BLAST databases from non-coral assemblies 
     base=$(basename ${infile} _contigs.fasta)
     # ITS2 genes
     # C. latusorum
     blastn -query ${prodir}/data/seqs/C_lat_ITS2_MW757029.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_lat_ITS2
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_lat_ITS2 | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_lat_ITS2.sorted
     # C. pacificum 
     blastn -query ${prodir}/data/seqs/C_pac_ITS2_MW757035.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_pac_ITS2
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_pac_ITS2 | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_C_pac_ITS2.sorted
     # D. glynni
     blastn -query ${prodir}/data/seqs/D_gly_ITS2_KY131788.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_D_gly_ITS2
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_D_gly_ITS2 | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/ITS2_blast_hits/${base}_D_gly_ITS2.sorted
done


# Put top hits from all samples into one text file
# python blast_sum.py