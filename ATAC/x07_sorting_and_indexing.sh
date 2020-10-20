############################################################################
### Sort and index bams
############################################################################

#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Sort all bams
for file in $(ls *deduplicated.bam)
do
  	base=$(basename $file "1_bismark_bt2_pe.deduplicated.bam")
    samtools sort -@ 6 -o ${base}deduplicated_sorted.bam ${file}
done

############################################################################

#!/bin/bash

#PBS -N indexing_bams
#PBS -l walltime=00:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.3.2

# Sort all bams
for file in $(ls *sorted.bam)
do
    samtools index ${file}
done




