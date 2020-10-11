#!/bin/bash

#PBS -N alignment
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16
#PBs -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

# Deal with the file names so they match
# Use basename to get the file name and then use a regex to make the symbols needed
# base=$(basename ${file})
# sample=${base:0:3} this gives the first three letters of the base
# replace this with regex to get the name between the hyphens
# this would also require rm the 'pristine' etc from the genome names


GENOME_PATH=/scratch/monoallelic/hm257/hm343_stuff/WGS/alternate_refs

for file in $(ls *R1.fq.gz)
do
	base=$(basename $file "R1.fq.gz")
    bowtie2 --very-sensitive -x ${GENOME_PATH}/${base} --threads 16 \
	-1 ${base}-trim_1.fq.gz \
	-2 ${base}-trim_2.fq.gz \
	-U ${base}-trim_unpaired_1.fq.gz \
	-U ${base}-trim_unpaired_2.fq.gz \
	-S ${base}.sam
done