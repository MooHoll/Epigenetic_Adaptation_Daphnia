    ## -------------------------------------------------------------------------
## Differential Methylation Between Temportal Populations: Daphnia.
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
library(methylKit)
library(readr)

## -------------------------------------------------------------------------
#!/bin/bash

#PBS -N making_methylkit_files
#PBS -l walltime=08:30:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

#cd $PBS_O_WORKDIR
#module load R/3.6.1
#R --save -q -f  making_files.R 

# Read in the .bam files and convert them to methylkit .txt files for easier future use
file.list <- list("eutrophic_12_3_sorted.bam","eutrophic_12_4_sorted.bam",
                  "eutrophic_12_5_1_sorted.bam","eutrophic_13_1_sorted.bam",
                  "eutrophic_13_2_sorted.bam","eutrophic_13_3_sorted.bam",
                  "eutrophic_13_5_1_sorted.bam","eutrophic_14_5_1_sorted.bam",
                  "eutrophic_15_5_1_sorted.bam","pesticide_6_2_sorted.bam",
                  "pesticide_6_3_sorted.bam","pesticide_7_3_sorted.bam",
                  "pesticide_7_5_sorted.bam","pesticide_7_5_4_sorted.bam",
                  "pesticide_8_5_3_sorted.bam","pesticide_9_20_sorted.bam",
                  "pesticide_9_5_1_sorted.bam","pesticide_9_5_3_sorted.bam",
                  "pesticide_9_6_sorted.bam","pristine_36_01_sorted.bam",
                  "pristine_36_02_sorted.bam","pristine_48_01_sorted.bam",
                  "pristine_48_02_sorted.bam","pristine_53_01_sorted.bam",
                  "pristine_54_01_sorted.bam","pristine_54_02_sorted.bam",
                  "pristine_74_01_sorted.bam","pristine_77_01_sorted.bam",
                  "pristine_88_01_sorted.bam","recovery_0_1_sorted.bam",
                  "recovery_0_2_sorted.bam","recovery_0_4_sorted.bam",
                  "recovery_1_2_sorted.bam","recovery_2_1_sorted.bam",
                  "recovery_2_5_11_sorted.bam","recovery_2_5_9_sorted.bam",
                  "recovery_3_5_1_sorted.bam","recovery_3_5_15_sorted.bam",
                  "recovery_3_5_2_sorted.bam","recovery_3_6_sorted.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                                               "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                                               "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                                               "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                                               "PP_LRV7_5_4", "PP_LRV7_4", "PP_LRV8_5_3",
                                               "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                                               "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                                               "PR_LR2_48_01", "PR_LR2_48_02", "PR_LR3_53_01",
                                               "PR_LR2_54_01", "PR_LR2_54_02", "PR_LR3_74_01",
                                               "PR_LR3_77_01", "PR_LR3_88_01", "CWP_LRV0_1",
                                               "CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                                               "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                                               "CWP_LRV3_5_15", "CWP_LRV3_5_1", "CWP_LRV3_5_2",
                                               "CWP_LRV3_6"),
                              treatment = c(rep(0,9), rep(1,10), rep(2,10), rep(3,11)),
                              assembly="daphmag_LRV0_1", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)


