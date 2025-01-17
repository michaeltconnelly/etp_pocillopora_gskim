#!/bin/bash
#
# Purpose: create concatenated symbiont genome reference for use with zooxtype.pl and ANGSD
# first, download original genome fasta and GFF files from source
# Cladocopium: http://symbs.reefgenomics.org/download/SymbC1.Genome.Scaffolds.fasta.gz, http://symbs.reefgenomics.org/download/SymbC1.Gene_Models.GFF3.gz
# Durusdinium: https://marinegenomics.oist.jp/symbd/download/102_symbd_genome_scaffold.fa.gz, https://marinegenomics.oist.jp/symbd/download/102_symbd.gff.gz
# rename to symC_genome and symD_genome

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
symdir="${mcs}/sequences/symbio"

# Create empty files to hold concatenated reference
touch ${symdir}/sym_catCD_genome.fasta ${symdir}/sym_catCD_genome.gff

# Append each reference genome into a concatenated reference
cat ${symdir}/cladocopium/symC1_genome.fasta ${symdir}/durusdinium/symD_genome.fasta >> ${symdir}/sym_catCD_genome.fasta

# Append each gff files into a concatenated reference
cat ${symdir}/cladocopium/symC1_genome.gff ${symdir}/durusdinium/symD_genome.gff >> ${symdir}/sym_catCD_genome.gff

# NOTES for use with zooxtype.pl
# -----------------------------------------------------------------------------------------------------
# Remove all but one fasta header for each reference genome, discard all newlines
# chr1, chr2
# C,    D,
grep \>.* ${symdir}/sym_catCD_genome.fasta

#>Smic.scaffold1|size3144590
#scaffold
#SymbC1
#sc0
# find and replace first instance with >chr1, etc.
sed '0,/\>Smic.*/{s/>Smic.*/\>chr1/}' ${symdir}/sym_cat_genome.fasta > ${symdir}/sym_cat_header1.fasta

# index refrences with BWA prior to using for alignment wrapper
#bwa index sym_catCD_genome.fasta
#bwa index symCDgenome.fasta