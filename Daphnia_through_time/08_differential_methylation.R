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

# Read in the .bam files and convert them to methylkit .txt files for easier future use
all_files <- read.csv()
file.list <- all_files

raw_data <- processBismarkAln(file.list,
                              sample.id = list(rep("Recovery", 10), rep("Pesticide", 10),
                                               rep("Eutrophic", 10), rep("Pristine", 10)),
                              assembly="ASM399081v1", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)

## -------------------------------------------------------------------------

# Read in combinations of files depending on which comparison is being made
file.list <- list("x_CpG.txt", "x_CpG.txt", "x_CpG.txt", "x_CpG.txt","x_CpG.txt",
                  "x_CpG.txt", "x_CpG.txt", "x_CpG.txt", "x_CpG.txt","x_CpG.txt",
                  "x_CpG.txt", "x_CpG.txt", "x_CpG.txt", "x_CpG.txt","x_CpG.txt",
                  "x_CpG.txt", "x_CpG.txt", "x_CpG.txt", "x_CpG.txt","x_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("j1nr","j1r","j5nr","j5r","j8nr","j8r"),
                     treatment = c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
                     assembly="ASM399081v1",
                     context="CpG")

## -------------------------------------------------------------------------

# Filter data for outliers and coverage
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)


## -------------------------------------------------------------------------

# Only text CpGs present in all samples (maybe too stingent, can adjust later)
meth_all_data <- unite(filtered_data, destrand=TRUE)

## -------------------------------------------------------------------------

# Filter sites using the MSC procedure so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)
df_meth_all$rownums <- row.names(df_meth_all)

source("/scratch/monoallelic/hm257/repro_methylation/bams/meth_calling_R_code/MSC.R")
source("/scratch/monoallelic/hm257/repro_methylation/bams/meth_calling_R_code/rateestimate.R")

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]
i <- df_meth_all[,26:27]
j <- df_meth_all[,29:30]
k <- df_meth_all[,32:33]
l <- df_meth_all[,35:36]
m <- df_meth_all[,38:39]
n <- df_meth_all[,41:42]
o <- df_meth_all[,44:45]
p <- df_meth_all[,47:48]
q <- df_meth_all[,50:51]
r <- df_meth_all[,53:54]
s <- df_meth_all[,56:57]
t <- df_meth_all[,59:60]

samples <- list(letters[1-20])

for(i in seq_along(samples)){
  colnames(samples[[i]]) <- c("CT","Ccount")
  MSCount <- MSC(samples[[i]], 1e-08)
  MSCresult <- MSCount$MSC
  pi <- MSCount$pi
  MSCrate <- rateestimate(MSCresult,pi)
  rate <- as.data.frame(MSCrate)
  MSCresult$row_nums <- row.names(MSCresult)
  label <- paste("MSCresult_meth_", i, sep="_")
  assign(label, subset(MSCresult, status == "methylated")) 
}

all_data <- rbind(paste0('MSCresult_meth__', 1:20))

meth_positions <- as.vector(as.numeric(unique(all_data$row_nums))) 

subset_methBase <- select (meth_all_data, meth_positions)

## -------------------------------------------------------------------------

# Save the dataframe for later use, including making nice plots
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="X_vs_Y_objectmethbase.txt", quote=F, row.names = F, sep = '\t')

## -------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("X_vs_Y_correlation.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("X_vs_Y_cluster.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("X_vs_Y_PCA.pdf")
PCASamples(subset_methBase)
dev.off()


## -------------------------------------------------------------------------

# Differential methylation
diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 3)
write.csv(diff_meth, file="X_vs_Y_all_tested_meth_sites_MSCfilter.csv")

diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.05)
write.csv(diff_meth_10, file="X_vs_Y_DMRs_min10percentDiff_qval0.05_MSCfilter.csv")


## -------------------------------------------------------------------------

# Annotating diff meth sites with gene IDs

genome_annotation<-read.csv.sql("file_with_genename_start_and_end.csv",
                                sql ="select * from file", sep=",",header = T)
colnames(diff_meth_10)[8]<-"meth_diff"


output <- sqldf("SELECT diff_meth_10.chr,
                diff_meth_10.start,
                diff_meth_10.end,
                diff_meth_10.pvalue,
                diff_meth_10.qvalue,
                diff_meth_10.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_meth_10 AS diff_meth_10
                LEFT JOIN genome_annotation AS ga 
                ON diff_meth_10.chr = ga.chr
                AND (diff_meth_10.start >= ga.start AND diff_meth_10.start <= ga.end)") 

output_dedup <- unique(output)

## -------------------------------------------------------------------------

# Write out all the info as well as 
write.table(output_dedup, file="X_vs_Y_diff_Meth_Genes_MSCfilter_with_geneID.csv",
            sep="\t", row.names=F, quote=F)

write.csv(output_dedup$geneID, file="X_vs_Y_list_DMGs_MSCfilter.csv")

