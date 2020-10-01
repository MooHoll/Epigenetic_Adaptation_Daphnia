#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=52:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bwa/0.7.17
module load samtools/1.9

# Define file paths
REF_FA=/scratch/monoallelic/hm343/daphnia/genome/PGA_assembly_DmagnaV3.fasta

# Remember to index the ref genome 1st: bwa index <genome_file>

# Align all fasta files to the reference
for file in $(ls *_1.fq.gz)
do
	base=$(basename ${file} "_1.fq.gz")
	bwa mem -t 16 ${REF_FA} ${base}_1.fq.gz ${base}_2.fq.gz > ${base}.sam
done
