#!/bin/bash
## from an input .bam file, randomly downsample to ds_ratio 
## of the original file
## Goal of this script is to reduce sequencing depth in 
## .pbcv1 .bam files to about 1000X, for computational reasons. 

## - fname_in = input filename
## - seed_ratio = seed for rng + downsampling ratio, separated by '.'

fname_in=$1
seed_ratio=$2

## load needed modules
module load gcc gdc samtools/1.3 

## create output filename
sname_out=$(echo ${fname_in} | sed 's/.bam/.ds1k.sam/g')
bname_out=$(echo ${fname_in} | sed 's/.bam/.ds1k.bam/g')

## copy header 
echo "###### DOWNSAMPLING ${fname_in} TO ${ds_ratio} OF ITS ORIGINAL SIZE ######"
samtools view -H ${fname_in} > ${sname_out}
## downsample reads
samtools view -s ${seed_ratio} ${fname_in} >> ${sname_out}
## transform to .bam
samtools view -b ${sname_out} > ${bname_out}

echo "###### DONE, CREATED ${bname_out} ######"
## remove .sam files
rm ${sname_out}