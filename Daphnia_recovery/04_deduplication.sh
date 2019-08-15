############################################################################
### Deduplication of .bam files
############################################################################

#!/bin/bash

#PBS -N deduplicating_bams
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.3.2
module load bismark/0.18.1

# Dedupliate all bam files, also gives a .txt report with how much data was removed
for file in $(ls *bam)
do
    deduplicate_bismark -p --bam ${file}
done
