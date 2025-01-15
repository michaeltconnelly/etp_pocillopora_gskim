#!/bin/bash

# symbiont barcode: psbA
# C. latusorum, C. pacificum, D.glynni

prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

for infile in $(ls ${prodir}/outputs/noncoral_spades/contigs/*_contigs.fasta)
do
     # Use previously made BLAST databases from non-coral assemblies 
     base=$(basename ${infile} _contigs.fasta)
     # psbA genes
     # C. latusorum
     blastn -query ${prodir}/data/seqs/C_lat_psbA_MW819757.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -perc_identity 97 -qcov_hsp_perc 50 -evalue 1e-6 -out ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_lat_psbA
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_lat_psbA | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_lat_psbA.sorted
     # C. pacificum 
     blastn -query ${prodir}/data/seqs/C_pac_psbA_MW861711.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -perc_identity 97 -qcov_hsp_perc 50 -evalue 1e-6 -out ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_pac_psbA
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_pac_psbA | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_C_pac_psbA.sorted
     # D. glynni
     blastn -query ${prodir}/data/seqs/D_gly_psbA_KY114947.fasta -db ${prodir}/outputs/noncoral_spades/blastdbs/${base} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -perc_identity 97 -qcov_hsp_perc 50 -evalue 1e-6 -out ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_D_gly_psbA
     sort -k3,3g -k12,12gr -k11,11g  ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_D_gly_psbA | sort -u -k1,1 --merge > ${prodir}/outputs/noncoral_spades/psbA_blast_hits/${base}_D_gly_psbA.sorted
done

# remove unsorted blast hits files
#rm  -- !(*.sorted)

# Put top hits from all samples into one text file
# python blast_sum_loci.py psbA