#!/bin/bash

#PBS -N calling_peaks
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load python/gcc/3.6.4

# needed to make a virtual environment:
# virtualenv --no-site-packages venv
# source /scratch/monoallelic/hm257/daphnia/bin/venv/bin/activate
# pip install numpy
# pip install MACS2

# to get chrom.sizes file
# pip install pyfaidx
# faidx PGA_assembly_DmagnaV3.fasta -i chromsizes > PGA_assembly_DmagnaV3.chrom.sizes
# somewhere along the line each chromosome gets renamed 1-71 urgh

bedGraphToBigWig=/scratch/monoallelic/hm257/daphnia/bin/bedGraphToBigWig
chrom_sizes=/scratch/monoallelic/hm257/daphnia/genome/PGA_assembly_DmagnaV3.chrom.sizes

source /scratch/monoallelic/hm257/daphnia/bin/venv/bin/activate

for file in $(ls *_sorted_filtered.bam)
do
  	base=$(basename ${file} "_sorted_filtered.bam")
    macs2 callpeak -t ${file} --format BAMPE --bdg --SPMR --gsize 12000000 --nolambda --keep-dup all --name ${base}.tmp.bg
    sort -k1,1 -k2,2n ${base}.tmp.bg_treat_pileup.bdg > ${base}.tmp.bedgraph
    ${bedGraphToBigWig} ${base}.tmp.bedgraph ${chrom_sizes} ${base}.bw
done
