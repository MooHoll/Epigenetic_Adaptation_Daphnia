## -------------------------------------------------------------------------
# Make one concensus peak file
## -------------------------------------------------------------------------
setwd("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/ATAC/2_peaks")

library(readr)
library(data.table)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*.narrowPeak")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
}

samples <- lapply(file.list, read_file1)

all_data <- rbindlist(samples)
all_data <- all_data[,c(1:3)] #71043

peaks <- all_data[!duplicated(all_data),] #58568 unique peaks
write.table(peaks, file ="concensus_peaks.bed", sep="\t", quote=F, col.names = F, row.names = F)
