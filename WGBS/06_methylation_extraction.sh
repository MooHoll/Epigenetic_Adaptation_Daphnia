############################################################################
### Extract methylation information to be used for weighted meth calculation
############################################################################

#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=10:00:00 + 10 hrs +
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run in current working directory
cd $PBS_O_WORKDIR

# Load modules
module load samtools/1.9

# Rename to match the new genome names
#for file in recovery_* ; do mv "${file}" "${file/recovery_/recovery-}"; done
#for file in pristine_* ; do mv "${file}" "${file/pristine_/pristine-}"; done
#for file in pesticide_* ; do mv "${file}" "${file/pesticide_/pesticide-}"; done
#for file in eutrophic_* ; do mv "${file}" "${file/eutrophic_/eutrophic-}"; done

for file in $(ls *.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    genome_name=$(echo ${base}| cut -d'-' -f 2)
    /scratch/monoallelic/hm257/hm343_stuff/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --comprehensive \
    --multicore 3 \
    --bedgraph \
    --cytosine_report \
    --genome_folder /scratch/monoallelic/hm257/hm343_stuff/WGS/alternate_refs/${genome_name} \
    ${file}
done


