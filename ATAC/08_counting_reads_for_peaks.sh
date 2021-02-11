#!/bin/bash

#PBS -N intersect_bed
#PBS -l walltime=03:00:00
#PBS -l vmem=100gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bedtools/2.28.0

# sort -k1,1 -k2,2n concensus_peaks.bed > concensus_peaks_sorted.bed 

for file in $(ls *bam)
do
	base=$(basename ${file} "_sorted_filtered.bam")
	bedtools intersect -wa -a concensus_peaks.bed -b ${file} -c > ${base}_peak_coverage.txt
done