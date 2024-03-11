#!/bin/bash
#./bash/ngsLD_scaffoldgroup_filter_prune_wrapper.sh
#purpose: create and submit jobs for ngsLD pruning on scaffold classes (eg. all scaffolds with 200-400 SNPs combined) and scaffolds >10K SNPs
#input: ANGSD SNP genotype likelihoods files, lists of sites across scaffolds 

#specify variable containing sequence file prefixes and directory paths
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"

# making a list of scaffold classes and scaffolds
groups=$(ls ${angsddir}/final_noclones/LD_pruning/sites | sed 's/_sites.txt//g')

#lets me know which files are being processed
echo "These are the groups:"
echo $groups

# loop to create and submit job scripts
for group in $groups
do \

# create empty ngsLD filtering job and pruning job files
echo "Creating jobs for ${group}"
#
NGSLD_JOB="${angsddir}/final_noclones/LD_pruning/jobs/ngsLD_${group}.job"
#
PRUNE_JOB="${angsddir}/final_noclones/LD_pruning/jobs/prune_${group}.job"

# populate ngsLD filtering job file
echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 8
#$ -q lThM.q
#$ -l mres=96G,h_data=12G,h_vmem=12G,himem
#$ -cwd
#$ -j y
#$ -N ngsLD_${group}_10k_filter
#$ -o ${angsddir}/final_noclones/LD_pruning/logs/ngsLD_${group}_10k_filter.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
 module load gcc/7.3/gsl
# ----------------Your Commands------------------- #" > $NGSLD_JOB

echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
NGSLD="/home/connellym/programs/ngsLD/ngsLD"
group=$1' >> $NGSLD_JOB

echo 'NS=$(wc -l ${angsddir}/final_noclones/LD_pruning/sites/${group}_sites.txt)
NB=$(wc -l ${angsddir}/final_noclones/final_noclones_bamfile.txt)' >> $NGSLD_JOB

echo '# run ngsLD command
$NGSLD --geno ${angsddir}/final_noclones/final_noclones_ibs05.geno.gz --probs 1 \
--n_ind $NB --n_sites $NS --pos ${angsddir}/final_noclones/LD_pruning/sites/${group}_sites.txt \
--max_kb_dist 10 \
--out ${angsddir}/final_noclones/LD_pruning/filter_out/${group}_LD.out --n_threads $NSLOTS --extend_out 1' >> $NGSLD_JOB

echo "# submit prune_graph job when finished
qsub $PRUNE_JOB $group" >> $NGSLD_JOB

echo 'echo = `date` job $JOB_NAME done' >> $NGSLD_JOB

# populate prune_graph job file

echo "#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q lThM.q
#$ -l mres=32G,h_data=8G,h_vmem=8G,himem
#$ -cwd
#$ -j y
#$ -N ngsLD_prune_${group}
#$ -o ${angsddir}/final_noclones/LD_pruning/logs/ngsLD_prune_${group}.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
# ----------------Your Commands------------------- #" > $PRUNE_JOB

echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim"
angsddir="/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/angsd"
PRUNE="/home/connellym/programs/prune_graph/target/release/prune_graph"
#
group=$1' >> $PRUNE_JOB

echo '# run prune_graph command
$PRUNE --in ${angsddir}/final_noclones/LD_pruning/filter_out/${group}_LD.out --header \
--weight-field "r2" --weight-filter "dist <= 10000 && r2 >= 0.5" \
--out ${angsddir}/final_noclones/prune_out/${group}_LD_prune.out -n $NSLOTS -v -v' >> $PRUNE_JOB

echo 'echo = `date` job $JOB_NAME done' >> $PRUNE_JOB

qsub $NGSLD_JOB $group #(calls prune job when finished)
#
done