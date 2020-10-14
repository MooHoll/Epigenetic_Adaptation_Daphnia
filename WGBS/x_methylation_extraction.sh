############################################################################
### Extract methylation information to be used for weighted meth calculation
############################################################################

#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=21:00:00
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
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    bismark_methylation_extractor -p \
    --comprehensive \
    --multicore 2 \
    --bedgraph \
    --cytosine_report \
    --genome_folder /scratch/monoallelic/hm257/daphnia/alternate_references/${base} \
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
module load R/3.4.1
library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", col_names=F)
}

samples <- lapply(file.list, read_file1)


for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,c(1,2,3,5)]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(i,"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}

# NOTE: will need to rename files with corresponding genotype, need to improve loop above for this.
#[1] "eutrophic_12_3_coverage.txt"   "eutrophic_12_4_coverage.txt"  
# [3] "eutrophic_12_5_1_coverage.txt" "eutrophic_13_1_coverage.txt"  
# [5] "eutrophic_13_2_coverage.txt"   "eutrophic_13_3_coverage.txt"  
# [7] "eutrophic_13_5_1_coverage.txt" "eutrophic_14_5_1_coverage.txt"
# [9] "eutrophic_15_5_1_coverage.txt" "pesticide_6_2_coverage.txt"   
#[11] "pesticide_6_3_coverage.txt"    "pesticide_7_3_coverage.txt"   
#[13] "pesticide_7_5_4_coverage.txt"  "pesticide_7_5_coverage.txt"   
#[15] "pesticide_8_5_3_coverage.txt"  "pesticide_9_20_coverage.txt"  
#[17] "pesticide_9_5_1_coverage.txt"  "pesticide_9_5_3_coverage.txt" 
#[19] "pesticide_9_6_coverage.txt"    "pristine_36_01_coverage.txt"  
#[21] "pristine_36_02_coverage.txt"   "pristine_48_01_coverage.txt"  
#[23] "pristine_48_02_coverage.txt"   "pristine_53_01_coverage.txt"  
#[25] "pristine_54_01_coverage.txt"   "pristine_54_02_coverage.txt"  
#[27] "pristine_74_01_coverage.txt"   "pristine_77_01_coverage.txt"  
#[29] "pristine_88_01_coverage.txt"   "recovery_0_1_coverage.txt"    
#[31] "recovery_0_2_coverage.txt"     "recovery_0_4_coverage.txt"    
#[33] "recovery_1_2_coverage.txt"     "recovery_2_1_coverage.txt"    
#[35] "recovery_2_5_11_coverage.txt"  "recovery_2_5_9_coverage.txt"  
#[37] "recovery_3_5_1_coverage.txt"   "recovery_3_5_15_coverage.txt" 
#[39] "recovery_3_5_2_coverage.txt"   "recovery_3_6_coverage.txt" 