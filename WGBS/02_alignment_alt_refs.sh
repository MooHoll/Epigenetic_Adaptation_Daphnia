############################################################################
### Alignment to alternate references
############################################################################

#!/bin/bash

#PBS -N genome_prep
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Prepare alt genomes
module load bowtie2/2.3.5.1
module load samtools/1.9

# Make a folder for each file and place it in
#for x in ./*.fasta; do
#  mkdir "${x%.*}" && cp "$x" "${x%.*}"
#done

# Run genome prep on folder containing .fa of genome (takes ~3 hours)
for folder in */
do
    /scratch/monoallelic/hm257/hm343_stuff/bin/Bismark-0.22.3/bismark_genome_preparation \
    ./${folder}
done

# -----------------------------------
# Alignment to parental genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=00:03:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16
#PBS -q devel

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.2.9
module load samtools/1.3.2

for file in $(ls *1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    /scratch/monoallelic/hm257/hm343_stuff/bin/Bismark-0.22.3/bismark \
    --multicore 8 -o alignment_to_alt_refs \
    /scratch/monoallelic/hm257/hm343_stuff/WGS/alternate_refs/${base} \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done