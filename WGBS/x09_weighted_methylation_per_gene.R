## -------------------------------------------------------------------------
## Weighted Methylation per Gene: Daphnia (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3523709/)
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time")
library(sqldf)
library(readr)
library(doBy)


## -------------------------------------------------------------------------

# Make the file which has the total CpGs per gene information from Bismark output (script 05)
cpgs <- read.delim("./Useful_Dmag_files/total_cpgs_in_genome.txt", header=F)

# File with gene start and end 
genes <- read.delim("./Useful_Dmag_files/genes_with_start_and_end.txt", header=F)
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

# Read in sample methylation count files (made in script 05)
# Take 1.5hrs on alice2
file.list = list.files(("./"),pattern="*final_coverage.txt")

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
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/Alignment_oldref_Genome/weighted_meth_per_gene")
file.list = list.files(("./"),pattern="*with_weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each population
eutrophic <- samples[1:9]
pesticide <- samples[10:19]
pristine <- samples[20:29]
recovery <- samples[30:40]

library(dplyr)
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

# Get a list of methylated genes which occur in at least one population

all_meth <- read_delim("Dmag_melted_weighted_meth_by_pop.txt", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


all_meth <- all_meth[,1:5]

# NOTE: this is an arbituary cut off of calling a feature methylated
# Here we say 15% of reads across a feature means it is a 'methylated' feature.
hist(all_meth$eutrophic[all_meth$eutrophic > 0.15])
hist(all_meth$pesticide[all_meth$pesticide > 0.15])
hist(all_meth$pristine[all_meth$pristine > 0.15])
hist(all_meth$recovery[all_meth$recovery > 0.15])

meth_genes <- all_meth[(all_meth$recovery > 0.15 | all_meth$pesticide > 0.15 | all_meth$eutrophic > 0.15
                        | all_meth$pristine > 0.15),]

meth_genes <- meth_genes[!is.na(meth_genes$geneID),]

genes_only <- as.data.frame(meth_genes$geneID)
write.table(genes_only, file = "Dmag_methylated_genes.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

all_meth_melted <- read_delim("Dmag_weighted_meth_by_pop.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
plot_data <- all_meth_melted[all_meth_melted$geneID %in% genes_only$`meth_genes$geneID`,]

plot_data$origin[plot_data$origin == "eutrophic"] <- "Eutrophic"
plot_data$origin[plot_data$origin == "pesticide"] <- "Pesticide"
plot_data$origin[plot_data$origin == "pristine"] <- "Pristine"
plot_data$origin[plot_data$origin == "recovery"] <- "Recovery"

plot_data$log <- (plot_data$weightedMeth.mean*plot_data$weightedMeth.mean)

library(ggforce)

ggplot(plot_data, aes(x=origin, y=weightedMeth.mean, fill=origin))+
  geom_boxplot()+
  theme_bw()+
  xlab("Population")+
  ylab("Mean Weighted Methylation Level per Gene")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title=element_text(size = 18))+
  guides(fill=FALSE)+
  scale_fill_manual(breaks = c("Pristine","Eutrophic","Pesticide","Recovery"),
                    values=c("#339966","#6600CC","#FFCC00", "#CC0066"))+
  scale_x_discrete(breaks = c("Pristine","Eutrophic","Pesticide","Recovery"),
                   limits = c("Pristine","Eutrophic","Pesticide","Recovery"))
  #facet_zoom(ylim = c(0.15, 0.3))
  


