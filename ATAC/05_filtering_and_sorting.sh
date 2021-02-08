#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Filter for mapped reads with a quality of at least 10 and sort and index the files
for file in $(ls *.bam)
do
  	base=$(basename ${file} ".bam")
    samtools view -b -F 4 -@ 8 ${file} | \
        samtools sort -@ 8 -m 3000000 -T ${base}.tmp.bam - | \
        samtools view -@ 8 -b -q 10 > ${base}_sorted_filtered.bam
    samtools index -@ 8 ${base}_sorted_filtered.bam
done