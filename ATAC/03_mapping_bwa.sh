#!/bin/bash

#PBS -N alignment
#PBS -l walltime=43:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load bwa/0.7.17
module load samtools/1.9 

# Index genome first
# bwa index PGA_assembly_DmagnaV3.fasta

genome=/scratch/monoallelic/hm257/daphnia/genome/PGA_assembly_DmagnaV3.fasta

for file in $(ls *1.fq.gz)
do
  	base=$(basename ${file} "_Rtrim_1.fq.gz")
    bwa mem -t 16 ${genome} ${base}_Rtrim_1.fq.gz ${base}_Rtrim_2.fq.gz | \
        samtools view -b -F 4 -@ 16 - | \
        samtools sort -@ 16 -m 3000000 -T ${base}.tmp.bam | \
        samtools view -@ 16 -b -q 10 > ${base}_sorted_filtered.bam
    samtools index -@ 16 ${base}_sorted_filtered.bam
done