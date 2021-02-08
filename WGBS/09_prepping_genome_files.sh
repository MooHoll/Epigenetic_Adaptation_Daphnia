############################################################################
### Prepare files
############################################################################

# Use one genome-wide cytosine report to create a file which consists of 
# just the scaffold name and the CpG position in a text file
cut -f1,2 eutrophic_12_3_1_bismark_bt2_pe.deduplicated.CpG_report.txt > total_cpgs_in_genome.txt


# From coverage files need files to look like:
# chr, position, total coverage, count cystosines
gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done

# R
module load R/3.6.1
library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read.delim(x, sep="\t", header=F)
}

samples <- lapply(file.list, read_file1)

sample.id <- list("EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                 "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                 "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                 "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                 "PP_LRV7_5_4", "PP_LRV7_5", "PP_LRV8_5_3",
                 "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                 "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                 "PR_LR2_48_01", "PR_LR2_48_02", "PR_LR3_53_01",
                 "PR_LR2_54_01", "PR_LR2_54_02", "PR_LR3_74_01",
                 "PR_LR3_77_01", "PR_LR3_88_01", "CWP_LRV0_1",
                 "CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                 "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                  "CWP_LRV3_5_1", "CWP_LRV3_5_15","CWP_LRV3_5_2",
                 "CWP_LRV3_6")
names(samples) <- sample.id


for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,c(1,2,3,5)]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(names(samples[i]),"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}