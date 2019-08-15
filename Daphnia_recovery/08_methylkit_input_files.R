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
#PBS -l walltime=05:20:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

#cd $PBS_O_WORKDIR
#module load R/3.5.1
#R --save -q -f  making_files.R 

# Read in the .bam files and convert them to methylkit .txt files for easier future use
file.list <- list("36_01_C4_L7_deduplicated_sorted.bam","36_01_C6_L7_deduplicated_sorted.bam","36_01_P4_L7_deduplicated_sorted.bam",
                  "36_01_P6_L8_deduplicated_sorted.bam","36_01_R6_L8_deduplicated_sorted.bam","53_01_C4_L8_deduplicated_sorted.bam",
                  "53_01_C6_L8_deduplicated_sorted.bam","53_01_P4_L8_deduplicated_sorted.bam","53_01_P6_L8_deduplicated_sorted.bam",
                  "53_01_R6_L8_deduplicated_sorted.bam","54_01_C4_L8_deduplicated_sorted.bam","54_01_C6_L8_deduplicated_sorted.bam",
                  "54_01_P4_L8_deduplicated_sorted.bam","54_01_P6_L8_deduplicated_sorted.bam","54_01_R6_L8_deduplicated_sorted.bam",
                  "6_2_C4_L6_deduplicated_sorted.bam","6_2_C6_L6_deduplicated_sorted.bam","6_2_P4_L6_deduplicated_sorted.bam",
                  "6_2_P6_L6_deduplicated_sorted.bam","6_2_R6_L6_deduplicated_sorted.bam","7_5_C4_L6_deduplicated_sorted.bam",
                  "7_5_C6_L7_deduplicated_sorted.bam","7_5_P4_L7_deduplicated_sorted.bam","7_5_P6_L7_deduplicated_sorted.bam",
                  "7_5_R6_L7_deduplicated_sorted.bam","9_20_C4_L7_deduplicated_sorted.bam","9_20_C6_L7_deduplicated_sorted.bam",
                  "9_20_P4_L7_deduplicated_sorted.bam","9_20_P6_L7_deduplicated_sorted.bam","9_20_R6_L7_deduplicated_sorted.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("36_01_C4","36_01_C6","36_01_P4",
                                            "36_01_P6","36_01_R6","53_01_C4",
                                            "53_01_C6","53_01_P4","53_01_P6",
                                            "53_01_R6","54_01_C4","54_01_C6",
                                            "54_01_P4","54_01_P6","54_01_R6",
                                            "6_2_C4","6_2_C6","6_2_P4",
                                            "6_2_P6","6_2_R6","7_5_C4",
                                            "7_5_C6","7_5_P4","7_5_P6",
                                            "7_5_R6","9_20_C4","9_20_C6",
                                            "9_20_P4","9_20_P6","9_20_R6"),
                              treatment = c(rep(1:5,6)),
                              assembly="daphmag_2.4", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)


