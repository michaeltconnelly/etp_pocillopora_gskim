#!/bin/bash
#./bash/symC_psbA_consensus_seqs.sh
#purpose: generate consensus Cladocopium psbA sequences

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# change to correct working directory to run for loops
cd ${prodir}/outputs/symC_psbA/

#   identify reference with highest number of reads mapped 
for i in *.filtered100-90.bam
do
    echo -ne "$i\t" ; samtools view $i | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 | head -1 | awk '{print $2"\t"$1}'
done > Summary_mapping_Cladocopium_psbA_100-90.tab

#   sort bam files 
cat Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c ; do samtools sort -o ${prodir}/outputs/symC_psbA/${a%.filtered100-90.bam}.sorted.bam $a; done

#   index sorted bam files 
for i in *.sorted.bam
do
    samtools index $i
done

#   Generate text pileup output
cat Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c
do
    samtools mpileup -Aa -f ${prodir}/data/seqs/cladocopium_psbA.fasta -r $b ${a%.filtered100-90.bam}.sorted.bam -o ${a%.filtered100-90.bam}.sorted.mpileup
done

# building consensus (minimum coverage of 2) 
for i in *.sorted.mpileup
do
    perl ${prodir}/bash/MpileuptoConsensus.pl -mpileup $i -Mincov 2 -out ${i}.consensus.fa
done

# keep consensus sequences in one fasta file 
for i in *sorted.mpileup.consensus.fa
do
    echo -e ">${i%_C_*}"; tail -n +2 $i
done > All_psbA-Cladocopium100-90_psbA.aln.mpileup.consensus.fa

