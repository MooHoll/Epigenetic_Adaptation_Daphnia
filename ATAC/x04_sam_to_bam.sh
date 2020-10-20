#!/bin/bash

#PBS -N sam_to_bam
#PBS -l walltime=12:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

for file in $(ls *.sam)
do
	base=$(basename $file ".sam")
    samtools view -S -b ${base}.sam > ${base}.bam
done

