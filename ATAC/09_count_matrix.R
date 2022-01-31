#------------------------------------------------------------------
# Differential open chromatin 
#------------------------------------------------------------------

# on ALICE
library(readr)
library(tidyr)

#------------------------------------------------------------------
# Make sample metadata
file_name <- list.files("./", pattern="*peak_coverage.txt")
sample_info <- data.frame(file_name)
sample_info$data <- sample_info$file_name
sample_info <- separate(sample_info, data, c("condition","genotype","bio_replicate","tech_replicate"),
                        sep = "-", fill = 'right')
sample_info$bio_replicate <- gsub("_peak_coverage.txt","",sample_info$bio_replicate)
sample_info$tech_replicate <- gsub("_peak_coverage.txt","",sample_info$tech_replicate)
sample_info$sample <- paste0(sample_info$genotype, "-", sample_info$bio_replicate, "-", sample_info$tech_replicate)
sample_info$sample <- gsub("-NA","",sample_info$sample)
rownames(sample_info) <- sample_info$sample

sample_info$population <- c(rep("CWP",9), rep("EP", 26), rep("CWP",20), rep("PR",7),
                            rep("CWP", 3), rep("PR",17), rep("PP",9), rep("PR",5),
                            rep("PP", 9), rep("PR",3),
                            rep("PP", 12), rep("PR",2))

write.table(sample_info, file="ATAC_sample_matrix.txt", sep='\t', quote = F, col.names = T, row.names = F)


#------------------------------------------------------------------
# Make count matrix

# Need to make a count matrix which consists of the row names being the peak IDs,
# each column name is the sample_id and each field is the count per sample per peak
files <- file.path(".", paste0(sample_info$file_name))
samples <- lapply(files, function(x)read.table(x, header=F))

names(samples) <- sample_info$sample

for (i in seq_along(samples)) { 
  samples[[i]]$peak <- paste0(samples[[i]]$V1,"-",samples[[i]]$V2, "-", samples[[i]]$V3)
  samples[[i]] <- as.data.frame(samples[[i]][,4:5])
  colnames(samples[[i]]) <- c(paste0(names(samples[i])),"peak")
}

all = as.data.frame(Reduce(function(...) merge(..., all=T, by = "peak"), samples))

write.table(all, file="ATAC_count_matrix.txt", sep='\t', quote = F, col.names = T, row.names = F)

