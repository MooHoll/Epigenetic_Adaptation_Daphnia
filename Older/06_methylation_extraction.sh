############################################################################
### Extract methylation information to be used for weighted meth calculation
############################################################################

#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=20:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run in current working directory
cd $PBS_O_WORKDIR

# Load modules
module load samtools/1.9


for file in $(ls *deduplicated.bam)
do
    base=$(basename ${file} "_bismark_bt2_pe.deduplicated.bam")
    /scratch/monoallelic/hm257/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --comprehensive \
    --multicore 3 \
    --bedgraph \
    --cytosine_report \
    --genome_folder /scratch/monoallelic/hm257/daphnia/genome \
    ${file}
done


