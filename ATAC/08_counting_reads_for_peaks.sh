#!/bin/bash

#PBS -N intersect_bed
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bedtools/2.28.0

for file in $(ls *bam)
do
	base=$(basename ${file} "_sorted_filtered.bam")
	bedtools intersect -wa -a concensus_peaks.bed -b ${file} -c -sorted > ${base}_peak_coverage.txt
done