#!/bin/bash

#PBS -N meth_extraction
#PBS -l walltime=36:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9

# make m-bias plots, reports and destranded file for methylkit
for file in $(ls pesticide_6_2*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_6_2 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_6_2 \
    ${file}
done

for file in $(ls pesticide_7_5*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_7_5 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_7_5 \
    ${file}
done

for file in $(ls pesticide_9_20*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_9_20 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pesticide_9_20 \
    ${file}
done

for file in $(ls pristine_36_01*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_36_01 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_36_01 \
    ${file}
done

for file in $(ls pristine_53_01*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_53_01 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_53_01 \
    ${file}
done

for file in $(ls pristine_54_01*deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_54_01 \
    ${file}

    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/daphnia/alternate_refs/pristine_54_01 \
    ${file}
done