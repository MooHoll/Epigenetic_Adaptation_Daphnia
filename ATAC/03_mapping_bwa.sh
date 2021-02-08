#!/bin/bash

#PBS -N alignment
#PBS -l walltime=36:00:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=24

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load bwa/0.7.17
module load samtools/1.9 

# Index genome first
# bwa index Daphnia_magna_LRV0_1.scaffolds.fa

genome=/scratch/monoallelic/hm257/daphnia/genome/Daphnia_magna_LRV0_1.scaffolds.fa

for file in $(ls *1.fq.gz)
do
  	base=$(basename ${file} "_Rtrim_1.fq.gz")
    bwa mem -t 24 -o ${base}.sam ${genome} ${base}_Rtrim_1.fq.gz ${base}_Rtrim_2.fq.gz 
done