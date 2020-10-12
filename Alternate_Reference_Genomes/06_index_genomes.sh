#!/bin/bash

#PBS -N index_genomes
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

for file in $(ls *.fasta)
do
    base=$(basename ${file} ".fasta")
    bowtie2-build --threads 16 ${file} ${base}
done

# Make a folder for each file and place it in
#for x in ./*.fasta; do
#  mkdir "${x%.*}" && mv "$x" "${x%.*}"
#done