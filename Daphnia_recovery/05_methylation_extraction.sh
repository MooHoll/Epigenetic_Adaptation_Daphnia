############################################################################
### Extract methylation information to be used for weighted meth calculation
############################################################################

#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=30:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run in current working directory
cd $PBS_O_WORKDIR

# Load modules
module load bismark/0.18.1
module load samtools/1.3.2

for file in $(ls *.bam)
do
    bismark_methylation_extractor -p \
    --comprehensive \
    --multicore 2 \
    --bedgraph \
    --cytosine_report \
    --genome_folder /scratch/monoallelic/hm257/daphnia/genome \
    ${file}
done


############################################################################
### Prepare files
############################################################################

# Use one genome-wide cytosine report to create a file which consists of 
# just the scaffold name and the CpG position in a text file
cut -f1,2 trim_0_1_L3_1_bismark_bt2_pe.deduplicated.CpG_report.txt > total_cpgs_in_genome.txt


# From coverage files need files to look like:
# chr, position, total coverage, count cystosines
gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done

# R
library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", col_names=F)
}

samples <- lapply(file.list, read_file1)

for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,-4]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(i,"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}

# NOTE: will need to rename files with corresponding sample name, should improve loop above for this.
# 36_01 C4
# 36_01 C6
# 36_01 P4
# 36_01 P6
# 36_01 R4
# 53_01
# 54_01
# 6_2
# 7_5
# 9_20









