############################################################################
### Download Daphnia genome and prepare
############################################################################

# Analysis note: phase reference genome with SNPs from individual genotypes when do real analysis!!!!!!

# One made in conjunction with Bham, beware 28000 scaffolds
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/632/505/GCA_001632505.1_daphmag2.4/GCA_001632505.1_daphmag2.4_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/632/505/GCA_001632505.1_daphmag2.4/GCA_001632505.1_daphmag2.4_genomic.gff.gz

# renamed .fna genome file to .fa and gunzip

module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2

# Run genome prep on folder containing .fa of genome
bismark_genome_prepration ./genome/


############################################################################
### Align trimmed samples to Daphnia genome
############################################################################

#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=32:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2

# Align all samples to the reference
REF_FA=/scratch/monoallelic/hm257/daphnia/new_genome

for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	bismark --multicore 3 -o alignment_daphnia ${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done

############################################################################
### Align trimmed samples to lambda genome
############################################################################

#!/bin/bash

#PBS -N alignment_to_lambda
#PBS -l walltime=04:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2

# Align all samples to the lambda
REF_FA=/scratch/monoallelic/hm257/lambda_genome

for file in $(ls *1.fq.gz)
do
  	base=$(basename $file "1.fq.gz")
    bismark --multicore 3 --prefix lambda \
    -o alignment_lambda ${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done