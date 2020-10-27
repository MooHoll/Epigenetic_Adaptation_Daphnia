############################################################################
### Sort and index bams
############################################################################

#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Sort all bams
for file in $(ls *.bam)
do
  	base=$(basename $file ".bam")
    samtools sort -@ 6 -o ${base}_sorted.bam ${file}
done

#Index all bams
for file in $(ls *sorted.bam)
do
    samtools index ${file}
done