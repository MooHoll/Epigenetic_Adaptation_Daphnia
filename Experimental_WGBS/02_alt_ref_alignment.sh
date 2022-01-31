#!/bin/bash

#PBS -N alignment_ER
#PBS -l walltime=30:00:00
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

# Genome prep
#for folder in */
#do
#    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
#done

# Align samples
for file in $(ls pesticide_6_2*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_6_2 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done

for file in $(ls pesticide_7_5*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_7_5 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done

for file in $(ls pesticide_9_20*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_9_20 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done

for file in $(ls pristine_36_01*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_36_01 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done

for file in $(ls pristine_53_01*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_53_01 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done

for file in $(ls pristine_54_01*1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
    --multicore 3 -o alignment_to_alt_refs \
    /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_54_01 \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done