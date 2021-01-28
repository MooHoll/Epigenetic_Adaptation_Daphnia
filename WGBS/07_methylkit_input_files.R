    ## -------------------------------------------------------------------------
## Differential Methylation Between Temportal Populations: Daphnia.
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

## -------------------------------------------------------------------------
#!/bin/bash

#PBS -N making_methylkit_files
#PBS -l walltime=08:30:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

#cd $PBS_O_WORKDIR
#module load R/3.5.1
#R --save -q -f  making_files.R 

# Read in the .bam files and convert them to methylkit .txt files for easier future use
file.list <- list("eutrophic_12_3_1_deduplicated_sorted.bam","eutrophic_12_4_1_deduplicated_sorted.bam",
                  "eutrophic_12_5_1_1_deduplicated_sorted.bam","eutrophic_13_1_1_deduplicated_sorted.bam",
                  "eutrophic_13_2_1_deduplicated_sorted.bam","eutrophic_13_3_1_deduplicated_sorted.bam",
                  "eutrophic_13_5_1_1_deduplicated_sorted.bam","eutrophic_14_5_1_1_deduplicated_sorted.bam",
                  "eutrophic_15_5_1_1_deduplicated_sorted.bam","pesticide_6_2_1_deduplicated_sorted.bam",
                  "pesticide_6_3_1_deduplicated_sorted.bam","pesticide_7_3_1_deduplicated_sorted.bam",
                  "pesticide_7_5_1_deduplicated_sorted.bam","pesticide_7_5_4_1_deduplicated_sorted.bam",
                  "pesticide_8_5_3_1_deduplicated_sorted.bam","pesticide_9_20_1_deduplicated_sorted.bam",
                  "pesticide_9_5_1_1_deduplicated_sorted.bam","pesticide_9_5_3_1_deduplicated_sorted.bam",
                  "pesticide_9_6_1_deduplicated_sorted.bam","pristine_36_01_1_deduplicated_sorted.bam",
                  "pristine_36_02_1_deduplicated_sorted.bam","pristine_48_01_1_deduplicated_sorted.bam",
                  "pristine_48_02_1_deduplicated_sorted.bam","pristine_53_01_1_deduplicated_sorted.bam",
                  "pristine_54_01_1_deduplicated_sorted.bam","pristine_54_02_1_deduplicated_sorted.bam",
                  "pristine_74_01_1_deduplicated_sorted.bam","pristine_77_01_1_deduplicated_sorted.bam",
                  "pristine_88_01_1_deduplicated_sorted.bam","recovery_0_1_1_deduplicated_sorted.bam",
                  "recovery_0_2_1_deduplicated_sorted.bam","recovery_0_4_1_deduplicated_sorted.bam",
                  "recovery_1_2_1_deduplicated_sorted.bam","recovery_2_1_deduplicated_sorted.bam",
                  "recovery_2_5_11_1_deduplicated_sorted.bam","recovery_2_5_9_1_deduplicated_sorted.bam",
                  "recovery_3_5_1_1_deduplicated_sorted.bam","recovery_3_5_15_1_deduplicated_sorted.bam",
                  "recovery_3_5_2_1_deduplicated_sorted.bam","recovery_3_6_1_deduplicated_sorted.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("eutrophic_12_3","eutrophic_12_4",
                                               "eutrophic_12_5_1","eutrophic_13_1",
                                               "eutrophic_13_2","eutrophic_13_3",
                                               "eutrophic_13_5_1","eutrophic_14_5_1",
                                               "eutrophic_15_5_1","pesticide_6_2",
                                               "pesticide_6_3","pesticide_7_3",
                                               "pesticide_7_5","pesticide_7_5_4",
                                               "pesticide_8_5_3","pesticide_9_20",
                                               "pesticide_9_5_1","pesticide_9_5_3",
                                               "pesticide_9_6","pristine_36_01",
                                               "pristine_36_02","pristine_48_01",
                                               "pristine_48_02","pristine_53_01",
                                               "pristine_54_01","pristine_54_02",
                                               "pristine_74_01","pristine_77_01",
                                               "pristine_88_01","recovery_0_1",
                                               "recovery_0_2","recovery_0_4",
                                               "recovery_1_2","recovery_2_1",
                                               "recovery_2_5_11","recovery_2_5_9",
                                               "recovery_3_5_1","recovery_3_5_15",
                                               "recovery_3_5_2","recovery_3_6"),
                              treatment = c(rep(1,9), rep(2,10), rep(3,10), rep(4,11)),
                              assembly="daphmag_2.4", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)


