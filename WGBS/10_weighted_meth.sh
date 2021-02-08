#!/bin/bash

#PBS -N weighted_meth
#PBS -l walltime=02:00:00
#PBS -l vmem=100gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=12

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load R/3.6.1 

# NOTE as there are so many files and they are huge, break this into chunks of 10 files at a time

R --save -q -f weighted_meth.R