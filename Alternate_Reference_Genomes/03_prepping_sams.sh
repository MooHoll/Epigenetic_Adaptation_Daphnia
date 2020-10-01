#-----------------------------------------------
# Getting files ready for SNP calling
#-----------------------------------------------

for file in $(ls *sam)
do
samtools flagstat ${file}
done

# -------------------------------------------

# In order to remove duplicates with the new samtools you need to sort the .bam by 
# name, run fixmate, then sort by coordinate, then mark and remove duplicates ... urgh

#!/bin/bash

#PBS -N sam_to_bam
#PBS -l walltime=10:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Make a bam and sort by name (pipes didn't work because of some weird thing with the sorting)
for file in $(ls *sam)
do
    base=$(basename $file ".sam")
    samtools sort -n -O bam -o ${base}_sorted.bam ${file}
done

# Run fixmate
for file in $(ls *bam)
do
    base=$(basename $file "_sorted.bam")
    samtools fixmate -m ${file} ${base}_fixname.bam
done

# Sort again by coordinate now
for file in $(ls *_fixname.bam)
do
    base=$(basename $file "_fixname.bam")
    samtools sort -o ${base}_sorted_coord.bam ${file}
done

# Mark and remove duplicates
for file in $(ls *_sorted_coord.bam)
do
    base=$(basename $file "_sorted_coord.bam")
    samtools markdup -r -s ${file} ${base}_nodups_final.bam
done

rm *fixmate.bam
rm *numbersort.bam
rm *sorted_coord.bam