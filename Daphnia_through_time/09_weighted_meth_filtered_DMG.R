## -------------------------------------------------------------------------
## Filter diff meth CpG genes for min weighted meth diff of 10%
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
library(readr)
library(reshape2)


## -------------------------------------------------------------------------

# Get percentage difference in methylation between population comparisons for each gene

weighted_meth<- read_delim("Dmag_weighted_meth_by_pop.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

recovery_info <- subset(weighted_meth, origin == "recovery")
colnames(recovery_info)[3] <- "recovery"
recovery_info <- recovery_info[,-2]

pesticide_info <- subset(weighted_meth, origin == "pesticide")
colnames(pesticide_info)[3] <- "pesticide"
pesticide_info <- pesticide_info[,-2]

eutrophic_info <- subset(weighted_meth, origin == "eutrophic")
colnames(eutrophic_info)[3] <- "eutrophic"
eutrophic_info <- eutrophic_info[,-2]

pristine_info <- subset(weighted_meth, origin == "pristine")
colnames(pristine_info)[3] <- "pristine"
pristine_info <- pristine_info[,-2]


Pris_Eutro <- merge(pristine_info, eutrophic_info, by = "geneID")
Pris_Eutro$diff <- 100*((Pris_Eutro$pristine - Pris_Eutro$eutrophic)/Pris_Eutro$eutrophic)

Eutro_Pest <- merge(eutrophic_info, pesticide_info, by = "geneID")
Eutro_Pest$diff <- 100*((Eutro_Pest$eutrophic - Eutro_Pest$pesticide)/Eutro_Pest$pesticide)

Pest_Recov <- merge(pesticide_info, recovery_info, by = "geneID")
Pest_Recov$diff <- 100*((Pest_Recov$pesticide - Pest_Recov$recovery)/Pest_Recov$recovery)


## -------------------------------------------------------------------------

# Filter for min 10% difference 

Pris_Eutro <- Pris_Eutro[!is.na(Pris_Eutro$diff),]
Pris_Eutro$diff[Pris_Eutro$diff==Inf] <- 100
Pris_Eutro_filtered_data <- Pris_Eutro[(Pris_Eutro$diff > 10 | Pris_Eutro$diff < -10),]

Eutro_Pest <- Eutro_Pest[!is.na(Eutro_Pest$diff),]
Eutro_Pest$diff[Eutro_Pest$diff==Inf] <- 100
Eutro_Pest_filtered_data <- Eutro_Pest[(Eutro_Pest$diff > 10 | Eutro_Pest$diff < -10),]

Pest_Recov <- Pest_Recov[!is.na(Pest_Recov$diff),]
Pest_Recov$diff[Pest_Recov$diff==Inf] <- 100
Pest_Recov_filtered_data <- Pest_Recov[(Pest_Recov$diff > 10 | Pest_Recov$diff < -10),]


## -------------------------------------------------------------------------

# Read in gene lists with a differentially methylated CpG
Pris_Eutro_CpG_genes <- read_csv("X_vs_Y_list_DMGs_MSCfilter.csv")
Eutro_Pest_CpG_genes <- read_csv("X_vs_Y_list_DMGs_MSCfilter.csv")
Pest_Recov_CpG_genes <- read_csv("X_vs_Y_list_DMGs_MSCfilter.csv")


## -------------------------------------------------------------------------

# Filter diff CpG genes for min diff of 10%
combined_Pris_Eutro <- Pris_Eutro_CpG_genes[(Pris_Eutro_CpG_genes$geneID %in% check$`unique(Pris_Eutro_filtered_data$geneID)`),]
combined_Pris_Eutro$hypermethylated <- "Eutrophic"
combined_Pris_Eutro$hypermethylated[combined_Pris_Eutro$diff > 0] <- "Pristine"
write.table(combined_Pris_Eutro, file="Pris_Eutro_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_Pris_Eutro$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="Pris_Eutro_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_Eutro_Pest <- Eutro_Pest_CpG_genes[(Eutro_Pest_CpG_genes$geneID %in% check$`unique(Eutro_Pest_filtered_data$geneID)`),]
combined_Eutro_Pest$hypermethylated <- "Pesticde"
combined_Eutro_Pest$hypermethylated[combined_Eutro_Pest$diff > 0] <- "Eutrophic"
write.table(combined_Eutro_Pest, file="Eutro_Pest_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_Eutro_Pest$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="Eutro_Pest_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )



combined_Pest_Recov <- Pest_Recov_CpG_genes[(Pest_Recov_CpG_genes$geneID %in% check$`unique(Pest_Recov_filtered_data$geneID)`),]
combined_Pest_Recov$hypermethylated <- "Recovery"
combined_Pest_Recov$hypermethylated[combined_Pest_Recov$diff > 0] <- "Pesticde"
write.table(combined_Pest_Recov, file="Pest_Recov_diff_meth_10%_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )

genes <- as.data.frame(combined_Pest_Recov$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="Pest_Recov_diff_meth_weighted_genes_only.txt", sep="\t", quote = F, col.names = T, row.names = F )


## -------------------------------------------------------------------------

# See if more hypermethylated in one pop compared to another

pristine_hyper <- subset(combined_Pris_Eutro, diff > 0) # x
eutrophic_hyper <- subset(combined_Pris_Eutro, diff < 0) # x

eutrophic_hyper <- subset(combined_Eutro_Pest, diff > 0) # x
pesticide_hyper <- subset(combined_Eutro_Pest, diff < 0) # x

pesticide_hyper <- subset(combined_Pest_Recov, diff > 0) # x
recovery_hyper <- subset(combined_Pest_Recov, diff < 0) # x

# Goodness of fit
observed = c(x, x)    # observed frequencies 
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

