#!/bin/bash

#PBS -N binomial
#PBS -l walltime=02:00:00
#PBS -l vmem=100gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load R/3.6.1 

R --save -q -f binomial_methylation.R