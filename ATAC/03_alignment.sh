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
module load bowtie2/2.3.5.1
module load samtools/1.9

# To get this working need to remove all pop names from genome files
#for file in pristine* ; do mv "${file}" "${file/pristine_/}"; done
#for file in eutrophic* ; do mv "${file}" "${file/eutrophic_/}"; done
#for file in recovery* ; do mv "${file}" "${file/recovery_/}"; done
#for file in pesticide* ; do mv "${file}" "${file/pesticide_/}"; done

GENOME_PATH=/scratch/monoallelic/hm257/hm343_stuff/WGS/alternate_refs

for file in $(ls *trim_1.fq.gz)
do
	base=$(basename ${file})
	genome_match=$(echo ${base}| cut -d'-' -f 1)
	sample=$(basename ${file} "_trim_1.fq.gz")
    bowtie2 --very-sensitive -x ${GENOME_PATH}/${genome_match} --threads 16 \
	-1 ${sample}_trim_1.fq.gz \
	-2 ${sample}_trim_2.fq.gz \
	-U ${sample}_trim_unpaired_1.fq.gz \
	-U ${sample}_trim_unpaired_2.fq.gz \
	-S ${sample}.sam
done