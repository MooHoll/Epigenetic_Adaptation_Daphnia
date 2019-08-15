## -------------------------------------------------------------------------
## Weighted Methylation per Gene: Daphnia
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetic_Recovery/Sequencing_Data")
library(sqldf)
library(readr)
library(doBy)
library(dplyr)


## -------------------------------------------------------------------------

# Making the file which has the total CpGs per gene information

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
# File with gene start and end 
genes <- read.delim("genes_with_start_and_end.txt", header=F)
colnames(genes) <- c("chr","start","end","geneID")

output <- sqldf("SELECT sg.V1,
                sg.V2,
                fg.chr,
                fg.start,
                fg.end,
                fg.geneID
                FROM cpgs AS sg
                LEFT JOIN genes AS fg 
                ON sg.V1 = fg.chr
                AND (sg.V2 >= fg.start AND sg.V2 <= fg.end)")
output <- output[!is.na(output$geneID),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ geneID+chr+start+end, data = output, FUN=sum)

write.table(final, file="Dmag_genes_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)


## -------------------------------------------------------------------------
# Take 1.5hrs on alice2
# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Read in gene with start/end and total CpGs per gene
genes_with_total_cpgs <- read_table2("Dmag_genes_with_total_cpgs.txt")
colnames(genes_with_total_cpgs)[5] <- "cpg_count"

# Calculate weighted meth for each gene for each sample

  for(i in seq_along(samples)){
    df <- samples[[i]]
    output <- sqldf("SELECT sg.chr,
                    sg.cpg,
                    sg.count_c,
                    sg.total_coverage,
                    fg.chr,
                    fg.start,
                    fg.end,
                    fg.geneID,
                    fg.cpg_count
                    FROM df AS sg
                    LEFT JOIN genes_with_total_cpgs AS fg 
                    ON sg.chr = fg.chr
                    AND (sg.cpg >= fg.start AND sg.cpg <= fg.end)")
    output <- output[!is.na(output$geneID),]
    output <- output[,-c(1,2,5,6,7)]
    check <- summaryBy(total_coverage + count_c ~ geneID + cpg_count, data=output, FUN=sum) 
    check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
    myfile <- file.path("./", paste0(i,"_","with_weighted_meth.txt"))
    write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
  }


## -------------------------------------------------------------------------

# Make one file covering all samples

file.list = list.files(("./"),pattern="*with_weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each experiment
control_4 <- samples[c(1,6,11,16,21,26)]
control_6 <- samples[c(2,7,12,17,22,27)]
pesticide_4 <- samples[c(3,8,13,18,23,28)]
pesticide_6 <- samples[c(4,9,14,19,24,29)]
recovery_6 <- samples[c(5,10,15,20,25,30)]

for(i in seq_along(control_4)){
  control_4[[i]]$origin <- "control_4"
  control_4[[i]] <- control_4[[i]][,c(1,5,6)]
}
control_4_all <- as.data.frame(bind_rows(control_4))
control_4_merged <- summaryBy(weightedMeth ~ geneID + origin, data = control_4_all, FUN=mean)

for(i in seq_along(control_6)){
  control_6[[i]]$origin <- "control_6"
  control_6[[i]] <- control_6[[i]][,c(1,5,6)]
}
control_6_all <- as.data.frame(bind_rows(control_6))
control_6_merged <- summaryBy(weightedMeth ~ geneID + origin, data = control_6_all, FUN=mean)

for(i in seq_along(pesticide_4)){
  pesticide_4[[i]]$origin <- "pesticide_4"
  pesticide_4[[i]] <- pesticide_4[[i]][,c(1,5,6)]
}
pesticide_4_all <- as.data.frame(bind_rows(pesticide_4))
pesticide_4_merged <- summaryBy(weightedMeth ~ geneID + origin, data = pesticide_4_all, FUN=mean)

for(i in seq_along(pesticide_6)){
  pesticide_6[[i]]$origin <- "pesticide_6"
  pesticide_6[[i]] <- pesticide_6[[i]][,c(1,5,6)]
}
pesticide_6_all <- as.data.frame(bind_rows(pesticide_6))
pesticide_6_merged <- summaryBy(weightedMeth ~ geneID + origin, data = pesticide_6_all, FUN=mean)

for(i in seq_along(recovery_6)){
  recovery_6[[i]]$origin <- "recovery_6"
  recovery_6[[i]] <- recovery_6[[i]][,c(1,5,6)]
}
recovery_6_all <- as.data.frame(bind_rows(recovery_6))
recovery_6_merged <- summaryBy(weightedMeth ~ geneID + origin, data = recovery_6_all, FUN=mean)


all_data <- rbind(control_4_merged,control_6_merged,pesticide_4_merged,pesticide_6_merged,recovery_6_merged )
write.table(all_data, file="Dmag_weighted_meth_by_experiment.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")

library(reshape2)
all_data_melted <- dcast(all_data, geneID ~ origin)
write.table(all_data_melted, file="Dmag_melted_weighted_meth_by_experiment.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")



## -------------------------------------------------------------------------

# Make another file with the same information but also broken down by pop and experiment
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetic_Recovery/Sequencing_Data/Weighted_meth")
file.list = list.files(("./"),pattern="*with_weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each experiment
control_4_pris <- samples[c(1,6,11)]
control_4_pest <- samples[c(16,21,26)]
control_6_pris <- samples[c(2,7,12)]
control_6_pest <- samples[c(17,22,27)]
pesticide_4_pris <- samples[c(3,8,13)]
pesticide_4_pest <- samples[c(18,23,28)]
pesticide_6_pris <- samples[c(4,9,14)]
pesticide_6_pest <- samples[c(19,24,29)]
recovery_6_pris <- samples[c(5,10,15)]
recovery_6_pest <- samples[c(20,25,30)]

# Replace with above list name and run for each
for(i in seq_along(recovery_6_pest)){
  recovery_6_pest[[i]]$origin <- "recovery_6_pest"
  recovery_6_pest[[i]] <- recovery_6_pest[[i]][,c(1,5,6)]
}
recovery_6_pest_all <- as.data.frame(bind_rows(recovery_6_pest))
recovery_6_pest_merged <- summaryBy(weightedMeth ~ geneID + origin, data = recovery_6_pest_all, FUN=mean)


all_data <- rbind(control_4_pris_merged,control_4_pest_merged,control_6_pris_merged,control_6_pest_merged,
                  pesticide_4_pris_merged,pesticide_4_pest_merged, pesticide_6_pris_merged,pesticide_6_pest_merged,
                  recovery_6_pris_merged,recovery_6_pest_merged)
write.table(all_data, file="Dmag_weighted_meth_by_experiment_and_pop.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")

library(reshape2)
all_data_melted <- dcast(all_data, geneID ~ origin)
write.table(all_data_melted, file="Dmag_melted_weighted_meth_by_experiment_and_pop.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")


## -------------------------------------------------------------------------

# Get a list of methylated genes which occur in st least one pop in one experiment

all_meth <- read_delim("../Useful_Dmag_files/Dmag_melted_weighted_meth_by_experiment_and_pop.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# Higher than error rate = 4373 genes
meth_genes <- all_meth[(all_meth$control_4_pris > 0.05 | all_meth$control_4_pest > 0.05 |
                          all_meth$control_6_pris > 0.05 | all_meth$control_6_pest > 0.05 | 
                          all_meth$pesticide_4_pris > 0.05 | all_meth$pesticide_4_pest > 0.05 | 
                          all_meth$pesticide_6_pris > 0.05 | all_meth$pesticide_6_pest > 0.05 | 
                          all_meth$recovery_6_pris > 0.05 | all_meth$recovery_6_pest > 0.05),]

meth_genes <- meth_genes[!is.na(meth_genes$geneID),] # 4037

genes_only <- as.data.frame(meth_genes$geneID)
write.table(genes_only, file = "Dmag_methylated_genes.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

write.table(meth_genes, file = "Dmag_methylated_genes_with_exp_and_pop_weighted_meth.txt", sep="\t",
            col.names = T, row.names = F, quote = F)




