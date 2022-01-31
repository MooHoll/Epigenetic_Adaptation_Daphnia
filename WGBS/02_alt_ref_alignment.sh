#!/bin/bash

#PBS -N alignment_ETT
#PBS -l walltime=45:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

# Make a folder for each alternate ref fasta and place it in
#for x in ./*.fasta; do
#  mkdir "${x%.*}" && mv "$x" "${x%.*}"
#done

# Genome prep (2hrs)
#for folder in */
#do
#    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
#done

# Align samples
for file in $(ls *1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/${base} \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done