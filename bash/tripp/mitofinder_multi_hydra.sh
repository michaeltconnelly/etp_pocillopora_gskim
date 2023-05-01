#!/bin/sh

path="$1"
reference="$2"
repath=`echo ${path} | cut -d "/" -f 1,2,3,4,5`
mkdir -p ${repath}/data/results/mitofinder

 for x in ${path}/*_R1_PE_trimmed.fastq.gz ; do 
    sample=`echo ${x} | cut -d "_" -f 1,2,3,4,5 | cut -d "/" -f 8`

    qsub -o ${repath}/jobs/logs/${sample}_mitofinder.log \
    -wd ${repath}/data/results/mitofinder \
    -N ${sample}_mitofinder \
    mitofinder_multi.job ${path} ${sample} ${reference}
done

# Mitofinder defaults to put the output in the same folder from which the 
# program was run, and there doesn't appear to be a way to change this. So, I
# added "-wd", which tells hydra to run the program in whatever directory is
# listed, in this case mitofinder/results/. Since this would also cause the log
# files to be placed here, I have also changed the log output "-o" to jobs/logs/

