#!/bin/bash
#./bash/mitofinder_assembly_wrapper.sh
#purpose: annotate manually circularized mitochondrial genome assemblies as species representatives

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"

# making a list of sample names
set=$1
#files=$(ls /scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/data/raw/)
#samples=$(echo "$files" | cut -d . -f 1 | sort -u)
samples=$(cat ${prodir}/data/${set}_samples.txt)

#lets me know which files are being processed
echo "These are the samples for mitofinder:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThM.q
#$ -l mres=192G,h_data=12G,h_vmem=12G,himem" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
echo "#$ -j y
#$ -N ${sample}_mitofinder_assembly
#$ -o ${prodir}/bash/jobs/${sample}_mitofinder_assembly.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
echo "module load bioinformatics/mitofinder" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job

# circularize mitofinder outputs for representative samples
sc="/scratch/nmnh_corals/connellym/programs/Simple-Circularise/simple_circularise.py"

# manually circularize assembly using simple_circularise.py -- account for multiple contigs!
echo "python $sc ${prodir}/outputs/mitofinder/${sample}/${sample}_MitoFinder_megahit_mitfi_Final_Results/${sample}_mtDNA_contig.fasta \
${prodir}/outputs/mitofinder/${sample}/${sample}_MitoFinder_megahit_mitfi_Final_Results/${sample}_mtDNA_contig_circular.fasta"  >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job

# sed rename
sed -i "s/>/>$sample/" ${prodir}/outputs/mitofinder/${sample}/${sample}_MitoFinder_megahit_mitfi_Final_Results/${sample}_mtDNA_contig_circular.fasta

#   input mitofinder assembly re-annotation command
# pocillopora.gb consists of Pocillopora-only reference mitogenomes
echo "mitofinder \
-j ${sample} \
-o 4 \
-r ${mcs}/sequences/pocillopora.gb \
-a ${prodir}/outputs/mitofinder/${sample}/${sample}_MitoFinder_megahit_mitfi_Final_Results/${sample}_mtDNA_contig_circular.fasta \
--allow-intron --intron-size 12000 --adjust-direction" >> "${prodir}"/bash/jobs/${sample}_mitofinder_assembly.job
# allow intron in ND5, perform tRNA annotation with mitfi

#
echo 'echo '${sample}' successfully processed' >> "${prodir}"/bash/jobs/${sample}_mitofinder_assembly.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
# submit job
qsub -wd ${prodir}/outputs/mitofinder ${prodir}/bash/jobs/${sample}_mitofinder_assembly.job
sleep 0.5
#
done