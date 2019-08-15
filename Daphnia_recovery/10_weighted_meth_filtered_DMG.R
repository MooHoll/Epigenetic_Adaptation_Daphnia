## -------------------------------------------------------------------------
## Filter diff meth CpG genes for min weighted meth diff of 10%
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetic_Recovery/Sequencing_Data")
library(readr)
library(reshape2)


## -------------------------------------------------------------------------

# Get percentage difference in methylation between population comparisons for each gene

weighted_meth<- read_delim("./Useful_Dmag_files/Dmag_weighted_meth_by_experiment.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

control_4_info <- subset(weighted_meth, origin == "control_4")
colnames(control_4_info)[3] <- "control_4"
control_4_info <- control_4_info[,-2]

control_6_info <- subset(weighted_meth, origin == "control_6")
colnames(control_6_info)[3] <- "control_6"
control_6_info <- control_6_info[,-2]

pesticide_4_info <- subset(weighted_meth, origin == "pesticide_4")
colnames(pesticide_4_info)[3] <- "pesticide_4"
pesticide_4_info <- pesticide_4_info[,-2]

pesticide_6_info <- subset(weighted_meth, origin == "pesticide_6")
colnames(pesticide_6_info)[3] <- "pesticide_6"
pesticide_6_info <- pesticide_6_info[,-2]

recovery_6_info <- subset(weighted_meth, origin == "recovery_6")
colnames(recovery_6_info)[3] <- "recovery_6"
recovery_6_info <- recovery_6_info[,-2]


C4_vs_C6 <- merge(control_4_info, control_6_info, by = "geneID")
C4_vs_C6$diff <- 100*((C4_vs_C6$control_4 - C4_vs_C6$control_6)/C4_vs_C6$control_6)

P4_vs_P6 <- merge(pesticide_4_info, pesticide_6_info, by = "geneID")
P4_vs_P6$diff <- 100*((P4_vs_P6$pesticide_4 - P4_vs_P6$pesticide_6)/P4_vs_P6$pesticide_6)

C4_vs_P4 <- merge(control_4_info, pesticide_4_info, by = "geneID")
C4_vs_P4$diff <- 100*((C4_vs_P4$control_4 - C4_vs_P4$pesticide_4)/C4_vs_P4$pesticide_4)

C6_vs_P6 <- merge(control_6_info, pesticide_6_info, by = "geneID")
C6_vs_P6$diff <- 100*((C6_vs_P6$control_6 - C6_vs_P6$pesticide_6)/C6_vs_P6$pesticide_6)

C6_vs_R6 <- merge(control_6_info, recovery_6_info, by = "geneID")
C6_vs_R6$diff <- 100*((C6_vs_R6$control_6 - C6_vs_R6$recovery_6)/C6_vs_R6$recovery_6)

P6_vs_R6 <- merge(pesticide_6_info, recovery_6_info, by = "geneID")
P6_vs_R6$diff <- 100*((P6_vs_R6$pesticide_6 - P6_vs_R6$recovery_6)/P6_vs_R6$recovery_6)


## -------------------------------------------------------------------------

# Filter for min 10% difference 

C4_vs_C6 <- C4_vs_C6[!is.na(C4_vs_C6$diff),]
C4_vs_C6$diff[C4_vs_C6$diff==Inf] <- 100
C4_vs_C6_filtered_data <- C4_vs_C6[(C4_vs_C6$diff > 10 | C4_vs_C6$diff < -10),]

P4_vs_P6 <- P4_vs_P6[!is.na(P4_vs_P6$diff),]
P4_vs_P6$diff[P4_vs_P6$diff==Inf] <- 100
P4_vs_P6_filtered_data <- P4_vs_P6[(P4_vs_P6$diff > 10 | P4_vs_P6$diff < -10),]

C4_vs_P4 <- C4_vs_P4[!is.na(C4_vs_P4$diff),]
C4_vs_P4$diff[C4_vs_P4$diff==Inf] <- 100
C4_vs_P4_filtered_data <- C4_vs_P4[(C4_vs_P4$diff > 10 | C4_vs_P4$diff < -10),]

C6_vs_P6 <- C6_vs_P6[!is.na(C6_vs_P6$diff),]
C6_vs_P6$diff[C6_vs_P6$diff==Inf] <- 100
C6_vs_P6_filtered_data <- C6_vs_P6[(C6_vs_P6$diff > 10 | C6_vs_P6$diff < -10),]

C6_vs_R6 <- C6_vs_R6[!is.na(C6_vs_R6$diff),]
C6_vs_R6$diff[C6_vs_R6$diff==Inf] <- 100
C6_vs_R6_filtered_data <- C6_vs_R6[(C6_vs_R6$diff > 10 | C6_vs_R6$diff < -10),]

P6_vs_R6 <- P6_vs_R6[!is.na(P6_vs_R6$diff),]
P6_vs_R6$diff[P6_vs_R6$diff==Inf] <- 100
P6_vs_R6_filtered_data <- P6_vs_R6[(P6_vs_R6$diff > 10 | P6_vs_R6$diff < -10),]

## -------------------------------------------------------------------------

# Read in gene lists with a differentially methylated CpG
C4_vs_C6_CpG_genes <- read.csv("./Differential_meth/control_vs_control/C4_vs_C6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")
P4_vs_P6_CpG_genes <- read.csv("./Differential_meth/pest_vs_pest/P4_vs_P6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")
C4_vs_P4_CpG_genes <- read.csv("./Differential_meth/control_vs_pest_gen4/C4_vs_P4_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")
C6_vs_P6_CpG_genes <- read.csv("./Differential_meth/control_vs_pest_gen6/C6_vs_P6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")
C6_vs_R6_CpG_genes <- read.csv("./Differential_meth/control_vs_recovery/C6_vs_R6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")
P6_vs_R6_CpG_genes <- read.csv("./Differential_meth/pest_vs_recovery/P6_vs_R6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
                               sep="\t")


## -------------------------------------------------------------------------

# Filter diff CpG genes for min diff of 10%
combined_C4_vs_C6 <- C4_vs_C6_filtered_data[(C4_vs_C6_filtered_data$geneID %in% unique(C4_vs_C6_CpG_genes$geneID)),]
nrow(combined_C4_vs_C6)
combined_C4_vs_C6$hypermethylated <- "Control_6"
combined_C4_vs_C6$hypermethylated[combined_C4_vs_C6$diff > 0] <- "Control_4"
write.table(combined_C4_vs_C6, file="C4_vs_C6_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_C4_vs_C6$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="C4_vs_C6_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_P4_vs_P6 <- P4_vs_P6_filtered_data[(P4_vs_P6_filtered_data$geneID %in% unique(P4_vs_P6_CpG_genes$geneID)),]
nrow(combined_P4_vs_P6)
combined_P4_vs_P6$hypermethylated <- "Pest_6"
combined_P4_vs_P6$hypermethylated[combined_P4_vs_P6$diff > 0] <- "Pest_4"
write.table(combined_P4_vs_P6, file="P4_vs_P6_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_P4_vs_P6$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="P4_vs_P6_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_C4_vs_P4 <- C4_vs_P4_filtered_data[(C4_vs_P4_filtered_data$geneID %in% unique(C4_vs_P4_CpG_genes$geneID)),]
nrow(combined_C4_vs_P4)
combined_C4_vs_P4$hypermethylated <- "Pest_4"
combined_C4_vs_P4$hypermethylated[combined_C4_vs_P4$diff > 0] <- "Control_4"
write.table(combined_C4_vs_P4, file="C4_vs_P4_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_C4_vs_P4$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="C4_vs_P4_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_C6_vs_P6 <- C6_vs_P6_filtered_data[(C6_vs_P6_filtered_data$geneID %in% unique(C6_vs_P6_CpG_genes$geneID)),]
nrow(combined_C6_vs_P6)
combined_C6_vs_P6$hypermethylated <- "Pest_6"
combined_C6_vs_P6$hypermethylated[combined_C6_vs_P6$diff > 0] <- "Control_6"
write.table(combined_C6_vs_P6, file="C6_vs_P6_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_C6_vs_P6$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="C6_vs_P6_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_C6_vs_R6 <- C6_vs_R6_filtered_data[(C6_vs_R6_filtered_data$geneID %in% unique(C6_vs_R6_CpG_genes$geneID)),]
nrow(combined_C6_vs_R6)
combined_C6_vs_R6$hypermethylated <- "Recovery_6"
combined_C6_vs_R6$hypermethylated[combined_C6_vs_R6$diff > 0] <- "Control_6"
write.table(combined_C6_vs_R6, file="C6_vs_R6_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_C6_vs_R6$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="C6_vs_R6_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_P6_vs_R6 <- P6_vs_R6_filtered_data[(P6_vs_R6_filtered_data$geneID %in% unique(P6_vs_R6_CpG_genes$geneID)),]
nrow(combined_P6_vs_R6)
combined_P6_vs_R6$hypermethylated <- "Recovery_6"
combined_P6_vs_R6$hypermethylated[combined_P6_vs_R6$diff > 0] <- "Pest_6"
write.table(combined_P6_vs_R6, file="P6_vs_R6_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_P6_vs_R6$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="P6_vs_R6_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )


## -------------------------------------------------------------------------

# See if more hypermethylated in one pop compared to another

control_4_hyper <- subset(combined_C4_vs_C6, diff > 0) # 3
nrow(control_4_hyper)
control_6_hyper <- subset(combined_C4_vs_C6, diff < 0) # 14
nrow(control_6_hyper)
# X-squared = 7.1176, df = 1, p-value = 0.007633

pest_4_hyper <- subset(combined_P4_vs_P6, diff > 0) # 6
nrow(pest_4_hyper)
pest_6_hyper <- subset(combined_P4_vs_P6, diff < 0) # 26
nrow(pest_6_hyper)
# X-squared = 12.5, df = 1, p-value = 0.000407

control_4_hyper <- subset(combined_C4_vs_P4, diff > 0) # 23
nrow(control_4_hyper)
pest_4_hyper <- subset(combined_C4_vs_P4, diff < 0) # 14
nrow(pest_4_hyper)
# X-squared = 2.1892, df = 1, p-value = 0.139

control_6_hyper <- subset(combined_C6_vs_P6, diff > 0) # 4
nrow(control_6_hyper)
pest_6_hyper <- subset(combined_C6_vs_P6, diff < 0) # 12
nrow(pest_6_hyper)
# X-squared = 4, df = 1, p-value = 0.0455

control_6_hyper <- subset(combined_C6_vs_R6, diff > 0) # 5
nrow(control_6_hyper)
recovery_6_hyper <- subset(combined_C6_vs_R6, diff < 0) # 17
nrow(recovery_6_hyper)
# X-squared = 6.5455, df = 1, p-value = 0.01052

pest_6_hyper <- subset(combined_P6_vs_R6, diff > 0) # 9
nrow(pest_6_hyper)
recovery_6_hyper <- subset(combined_P6_vs_R6, diff < 0) # 10
nrow(recovery_6_hyper)
# X-squared = 0.052632, df = 1, p-value = 0.8185

# Goodness of fit
observed = c(9, 10)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 


## Maybe useful for the future
#male_hyper <- as.data.frame(male_hyper$geneID)
#colnames(male_hyper) <- "geneID"
#queen_hyper <- as.data.frame(queen_hyper$geneID)
#colnames(queen_hyper) <- "geneID"

#write.table(male_hyper, file="male_hyper_genes.txt", sep="\t", quote = F, col.names = T, row.names = F )
#write.table(queen_hyper, file="queen_hyper_genes.txt", sep="\t", quote = F, col.names = T, row.names = F )

