############################################################################
### Alignment to alternate references
############################################################################

#!/bin/bash

#PBS -N alignment_meth
#PBS -l walltime=36:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.2.9
module load samtools/1.3.2

# Need to rub bismark_prepare_genome command on the genome beforehand
# /scratch/monoallelic/hm257/bin/Bismark-0.22.3/bismark_genome_preparation ./genome

for file in $(ls *1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hm257/bin/Bismark-0.22.3/bismark \
    --multicore 8 -o alignment \
    /scratch/monoallelic/hm257/daphnia/genome \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done