#!/bin/bash

#PBS -N SNP_calling
#PBS -l walltime=30:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in (6hrs)
cd $PBS_O_WORKDIR 

# Load software needed
module load freebayes/1.1.0

# Define file paths
REF_FILE=/scratch/monoallelic/hm257/hm343_stuff/genome/PGA_assembly_DmagnaV3.fasta                                                                                      

# create the directory where the output files are to be written 
OUTPUT=vcf_files                                                                                                                                      
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

                                                                                                         
# Run freebayes for samples, min count 2 of alternative alleles, min 8 reads per SNP, ignore complex events, indels and mnps
# Takes ~ 1hr per sample urgh, can't parallelise easily
for file in $(ls *bam)
do
    freebayes \
        -f ${REF_FILE} \
        -C 2 \
        -! 8 \
        -u \
        -i \
        -X \
        -b ${file}\
        > ${OUTPUT}/${file}.freebayes.vcf
done

#-------------------------------------------------------------------------------

#!/bin/bash

#PBS -N snp_filtering
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR 

module load gcc/6.3.0
module load vcftools/0.1.14

# Takes 1min per sample (can run on login nodes)
for file in $(ls *.vcf)
do
    base=$(basename $file ".vcf")
    vcftools \
    --vcf ${file} \
    --minQ 20 \
    --min-meanDP 5 \
    --recode \
    --recode-INFO-all \
    --out ${base}
done