#!/bin/sh


path="$1"
mkdir ${path}/fastqc_analyses

 for x in ${path}/*.fastq.gz ; do 
  base=`echo ${x} | cut -d "_" -f 1,2,3,4,5 | cut -d "/" -f 8`
  

  qsub -o logs/${base}_fastqc.log \
  -N ${base}_fastqc \
  fastqc_multi.job ${path} ${base}

done

# This script runs a fastqc analysis on multiple fastq sequences, usually the
# output of a miseq/hiseq/novaseq run. It assumes the file names are the typical
# demultiplexed sample names downloaded by basespace. All files should be in the
# same folder usually PROJECT/data/raw. This file should be in PROJECT/jobs.
# For each sample file, it will output an html file that you can open with your
# browser and view visual represntations of your sample quality.
# It will also output a .zip file that contains text summaries of read quality.  
