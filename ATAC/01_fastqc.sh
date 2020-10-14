#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=08:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load fastqc/0.11.5

# create the directory where the output files are to be written 
OUTPUT=fastqc
if [ ! -d "$OUTPUT" ]; then
    mkdir -p ${OUTPUT}
fi

# Create a list of the files to be called
for file in $(ls *.fq.gz)
do
  	fastqc -t 14 -o ${OUTPUT} ${file}
done