## -------------------------------------------------------------------------
## Weighted Methylation per Gene: Daphnia
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
library(sqldf)
library(readr)
library(doBy)


## -------------------------------------------------------------------------

# Making the file which has the total CpGs per gene information

cpgs <- read.delim("destranded_cpg_positions.txt", header=F)
# File with gene start and end 
genes <- read.delim("x", header=T)

output <- sqldf("SELECT sg.V1,
                sg.V2,
                fg.chr,
                fg.start,
                fg.end,
                fg.feature,
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

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Read in gene with start/end and total CpGs per gene
genes_with_total_cpgs <- read_table2("Dmag_genes_with_total_cpgs.txt")
colnames(genes_with_total_cpgs)[5] <- "cpg_count"


## -------------------------------------------------------------------------

# Calculate weighted meth for each gene for each sample

for(i in seq_along(samples)){
  df <- samples[[i]]
  output <- sqldf("SELECT sg.chr,
                  sg.base,
                  sg.coverage,
                  sg.count_C,
                  fg.chr,
                  fg.start,
                  fg.end,
                  fg.geneID,
                  fg.cpg_count
                  FROM df AS sg
                  LEFT JOIN genes_with_total_cpgs AS fg 
                  ON sg.chr = fg.chr
                  AND (sg.base >= fg.start AND sg.base <= fg.end)")
  output <- output[!is.na(output$geneID),]
  output <- output[,-c(1,2,5,6,7)]
  check <- summaryBy(coverage + count_C ~ geneID + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_C.sum)/(check$cpg_count*check$coverage.sum)
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


# Make one dataframe for each population
recovery <- samples[1:10]
pesticide <- samples[11:20]
eutrophic <- samples[21:30]
pristine <- samples[31:40]

for(i in seq_along(recovery)){
  recovery[[i]]$origin <- "recovery"
  recovery[[i]] <- recovery[[i]][,c(1,5,6)]
}
recovery_all <- as.data.frame(bind_rows(recovery))
recovery_merged <- summaryBy(weightedMeth ~ geneID + origin, data = recovery_all, FUN=mean)


for(i in seq_along(pesticide)){
  pesticide[[i]]$origin <- "pesticide"
  pesticide[[i]] <- pesticide[[i]][,c(1,5,6)]
}
pesticide_all <- as.data.frame(bind_rows(pesticide))
pesticide_merged <- summaryBy(weightedMeth ~ geneID + origin, data = pesticide_all, FUN=mean)


for(i in seq_along(eutrophic)){
  eutrophic[[i]]$origin <- "eutrophic"
  eutrophic[[i]] <- eutrophic[[i]][,c(1,5,6)]
}
eutrophic_all <- as.data.frame(bind_rows(eutrophic))
eutrophic_merged <- summaryBy(weightedMeth ~ geneID + origin, data = eutrophic_all, FUN=mean)


for(i in seq_along(pristine)){
  pristine[[i]]$origin <- "pristine"
  pristine[[i]] <- pristine[[i]][,c(1,5,6)]
}
pristine_all <- as.data.frame(bind_rows(pristine))
pristine_merged <- summaryBy(weightedMeth ~ geneID + origin, data = pristine_all, FUN=mean)


all_data <- rbind(recovery_merged, pesticide_merged, eutrophic_merged, pristine_merged)
write.table(all_data, file="Dmag_weighted_meth_by_pop.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")


all_data_melted <- dcast(all_data, geneID ~ origin)
write.table(all_data_melted, file="Dmag_melted_weighted_meth_by_pop.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")


## -------------------------------------------------------------------------

# Get a list of methylated genes which occur in st least one population

all_meth <- read_delim("Dmag_melted_weighted_meth_by_pop.txt", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


all_meth <- all_meth[,1:5]
meth_genes <- all_meth[(all_meth$recovery > 0 | all_meth$pesticide > 0 | all_meth$eutrophic > 0
                        | all_meth$pristine > 0),]

meth_genes <- meth_genes[!is.na(meth_genes$geneID),]

genes_only <- as.data.frame(meth_genes$geneID)
write.table(genes_only, file = "Dmag_methylated_genes.txt", sep="\t",
            col.names = T, row.names = F, quote = F)



