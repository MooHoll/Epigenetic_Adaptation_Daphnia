############################################################################
### Fastqc and Trimming
############################################################################

#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load fastqc/0.11.5

# create the directory where the output files are to be written 
OUTPUT=fastqc
if [ ! -d "$OUTPUT" ]; then
    mkdir -p ${OUTPUT}
fi

# Create a list of the files to be called
for file in $(ls *.fq.gz)
do
  	fastqc -o $OUTPUT ${file}
done

############################################################################

#!/bin/bash

#PBS -N trimming
#PBS -l walltime=07:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed (pigz needed for parallel zipping)
module load cutadapt/1.11
module load pigz/2.3.3

###-------------------------------------
### For just the pristine samples:
###-------------------------------------

# Trim adapters (must be done before quality trimming)
for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	cutadapt \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o trim_${base}1.fq.gz \
	-p trim_${base}2.fq.gz \
	${base}1.fq.gz \
	${base}2.fq.gz
done

# Trim 10 bases from reads
for file in $(ls trim*1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	cutadapt -u 10 -u -10 -U 10 -U -10\
	-o trim2_${base}1.fq.gz \
	-p trim2_${base}2.fq.gz \
	${base}1.fq.gz \
	${base}2.fq.gz
done

############################################################################

#!/bin/bash

#PBS -N trimming
#PBS -l walltime=06:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed (pigz needed for parallel zipping)
module load cutadapt/1.11
module load pigz/2.3.3

###-------------------------------------
### For all other samples:
###-------------------------------------

for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	cutadapt -u 10 -u -10 -U 10 -U -10 \
	-o trim_${base}1.fq.gz \
	-p trim_${base}2.fq.gz \
	${base}1.fq.gz \
	${base}2.fq.gz
done



############################################################################

#!/bin/bash

#PBS -N fastqc_trimmed
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load fastqc/0.11.5

# create the directory where the output files are to be written 
OUTPUT=fastqc_trimmed
if [ ! -d "$OUTPUT" ]; then
    mkdir -p ${OUTPUT}
fi

# Create a list of the files to be called (time above for 24 files) 
for file in $(ls *.fq.gz)
do
  	fastqc -o $OUTPUT ${file}
done