## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(reshape2)

## -------------------------------------------------------------------------
# Read in weighted methylation files
file.list = list.files(("."),pattern="*weighted_meth.txt")

read_file1 <- function(x){
  read.delim(x, sep="\t", header = T)
}

samples <- lapply(file.list, read_file1)


sample.id = list("CWP_LRV0_1","CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                 "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                 "CWP_LRV3_5_15", "CWP_LRV3_5_1", "CWP_LRV3_5_2",
                 "CWP_LRV3_6","EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                 "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                 "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                 "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                 "PP_LRV7_5_4", "PP_LRV7_5", "PP_LRV8_5_3",
                 "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                 "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                 "PR_LR2_48_01", "PR_LR2_48_02",
                 "PR_LR2_54_01", "PR_LR2_54_02",  "PR_LR3_53_01",
                 "PR_LR3_74_01",
                 "PR_LR3_77_01", "PR_LR3_88_01")
names(samples) <- sample.id

# -----------------------------------------------------------------
# Make a dataframe which contains the information for each individual (this will be huge)
for (i in seq_along(samples)) { 
  samples[[i]] <- as.data.frame(samples[[i]][,-c(6:7)])
  colnames(samples[[i]]) <- c("feature","id","start","end","cpg_count",paste0(names(samples[i])))
}

all = as.data.frame(Reduce(function(...) merge(..., all=T, by = c("feature","id","start","end","cpg_count")), samples))

write.table(all, file="Dmagna_weighted_meth_all_samples.txt", sep='\t', quote = F, col.names = T, row.names = F)

# -----------------------------------------------------------------
# Make a dataframe which has the average weighted methylation level by population
# Read in weighted methylation files
file.list = list.files(("."),pattern="*weighted_meth.txt")

read_file1 <- function(x){
  read.delim(x, sep="\t", header = T)
}

samples <- lapply(file.list, read_file1)


sample.id = list("CWP_LRV0_1","CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                 "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                 "CWP_LRV3_5_15", "CWP_LRV3_5_1", "CWP_LRV3_5_2",
                 "CWP_LRV3_6","EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                 "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                 "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                 "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                 "PP_LRV7_5_4", "PP_LRV7_5", "PP_LRV8_5_3",
                 "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                 "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                 "PR_LR2_48_01", "PR_LR2_48_02",
                 "PR_LR2_54_01", "PR_LR2_54_02",  "PR_LR3_53_01",
                 "PR_LR3_74_01",
                 "PR_LR3_77_01", "PR_LR3_88_01")
names(samples) <- sample.id

CWP <- samples[1:11]
EP <- samples[12:20]
PP <- samples[21:30]
PR <- samples[31:40]

for(i in seq_along(CWP)){
  CWP[[i]]$population <- "CWP"
}
CWP_all <- as.data.frame(bind_rows(CWP))
CWP_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                          cpg_count + population, data = CWP_all, FUN=mean)

for(i in seq_along(EP)){
  EP[[i]]$population <- "EP"
}
EP_all <- as.data.frame(bind_rows(EP))
EP_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                              cpg_count + population, data = EP_all, FUN=mean)

for(i in seq_along(PP)){
  PP[[i]]$population <- "PP"
}
PP_all <- as.data.frame(bind_rows(PP))
PP_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                          cpg_count + population, data = PP_all, FUN=mean)

for(i in seq_along(PR)){
  PR[[i]]$population <- "PR"
}
PR_all <- as.data.frame(bind_rows(PR))
PR_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                          cpg_count + population, data = PR_all, FUN=mean)


all_data <- rbind(CWP_merged, EP_merged, PP_merged, PR_merged)

all_data2 <- dcast(all_data, feature + gene_id + start + end +
                     cpg_count ~ population, value.var = "weightedMeth.mean")

write.table(all_data2, file="Dmagna_weighted_meth_annotation_by_population.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")

# -----------------------------------------------------------------
#Add in the exon numbers using the sql trick

