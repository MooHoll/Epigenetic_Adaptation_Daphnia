## -------------------------------------------------------------------------
## Differential Methylation Between Control Gen 4 and Gen 6
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

## -------------------------------------------------------------------------

# Read in combinations of files depending on which comparison is being made
# It sucks but you have to name them rather than supplying a list with *

file.list <- list("36_01_C4_CpG.txt","36_01_C6_CpG.txt",
                  "53_01_C4_CpG.txt",
                  "53_01_C6_CpG.txt",
                  "54_01_C4_CpG.txt","54_01_C6_CpG.txt",
                  "6_2_C4_CpG.txt","6_2_C6_CpG.txt",
                  "7_5_C4_CpG.txt",
                  "7_5_C6_CpG.txt",
                  "9_20_C4_CpG.txt","9_20_C6_CpG.txt")


raw_data <- methRead(file.list,
                     sample.id = list("36_01_C4","36_01_C6",
                                      "53_01_C4",
                                      "53_01_C6",
                                      "54_01_C4","54_01_C6",
                                      "6_2_C4","6_2_C6",
                                      "7_5_C4",
                                      "7_5_C6",
                                      "9_20_C4","9_20_C6"),
                     treatment = c(rep(0:1, 6)),
                     assembly="daphmag_2.4",
                     context="CpG")

## -------------------------------------------------------------------------

# Filter data for outliers and coverage and normalise counts
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

normalized <- normalizeCoverage(filtered_data)

## -------------------------------------------------------------------------

# Only text CpGs present in all samples (maybe too stingent, can adjust later)
meth_all_data <- unite(normalized, destrand=TRUE)

## -------------------------------------------------------------------------

# Filter sites using the MSC procedure so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)
df_meth_all$rownums <- row.names(df_meth_all)

source("/scratch/monoallelic/hm257/MSC_scripts/MSC.R")
source("/scratch/monoallelic/hm257/MSC_scripts/rateestimate.R")

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]
j <- df_meth_all[,29:30]
k <- df_meth_all[,32:33]
l <- df_meth_all[,35:36]
m <- df_meth_all[,38:39]


samples <- list(a,b,c,d,e,f,g,h,j,k,l,m)

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

all_data <- rbind(MSCresult_meth__1,MSCresult_meth__2,MSCresult_meth__3,
                  MSCresult_meth__4,MSCresult_meth__5,MSCresult_meth__6,
                  MSCresult_meth__7,MSCresult_meth__8,MSCresult_meth__9,
                  MSCresult_meth__10,MSCresult_meth__11,MSCresult_meth__12)

meth_positions <- as.vector(as.numeric(unique(all_data$row_nums))) 

subset_methBase <- select (meth_all_data, meth_positions)

## -------------------------------------------------------------------------

# Save the dataframe for later use, including making nice plots
methBase_ob <- getData(subset_methBase) 
nrow(methBase_ob) #3938 Cpgs left for testing differential meth
write.table(methBase_ob, file="C4_vs_C6_objectmethbase.txt", quote=F, row.names = F, sep = '\t')

## -------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("C4_vs_C6_correlation.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("C4_vs_C6_cluster.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("C4_vs_C6_PCA.pdf")
PCASamples(subset_methBase)
dev.off()


## -------------------------------------------------------------------------

# Differential methylation
covariates <- data.frame(population = c("pristine","pristine","pristine","pristine","pristine","pristine",
                                        "pesticide","pesticide","pesticide","pesticide","pesticide","pesticide"),
                         genotype = c ("36_01", "36_01", "53_01","53_01","54_01", "54_01","6_2","6_2","7_5","7_5",
                                       "9_20","9_20"))

# Calculated as C6 - C4
diff_meth <- calculateDiffMeth(subset_methBase, covariates = covariates)
write.csv(diff_meth, file="C4_vs_C6_all_tested_meth_sites_MSCfilter.csv")

diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.05)
nrow(diff_meth_10) #number of differentially methylated CpGs
write.csv(diff_meth_10, file="C4_vs_C6_DMRs_min10percentDiff_qval0.05_MSCfilter.csv")


## -------------------------------------------------------------------------

# Annotating diff meth sites with gene IDs

genome_annotation<-read.csv.sql("genes_with_start_and_end.txt",
                                sql ="select * from file", sep="\t",header = F)
colnames(genome_annotation) <- c("chr","start","end","geneID")

diff_meth_10 <- getData(diff_meth_10)
colnames(diff_meth_10)[7]<-"meth_diff"


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

# Write out all the info as well as just a file with the gene names
write.table(output_dedup, file="C4_vs_C6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
            sep="\t", row.names=F, quote=F)

write.csv(output_dedup$geneID, file="C4_vs_C6_list_DMGs_MSCfilter.csv")

## -------------------------------------------------------------------------

# Forgot to remove rows with NA to see how many CpGs in gene and not
library(readr)
data1 <- read.csv("C4_vs_C6_diff_Meth_Genes_MSCfilter_with_geneID.csv", sep="\t")
data2 <- data1[!is.na(data1$geneID),]
nrow(data1)
nrow(data2)
write.table(data2, file="C4_vs_C6_diff_Meth_Genes_MSCfilter_with_geneID.csv",
            sep="\t", row.names=F, quote=F)

# 38 CpGs in genes and 7 not in genes



