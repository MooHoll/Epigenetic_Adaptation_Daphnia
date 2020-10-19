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
# Put in the correct comparison below, see script: 08_comparisons_diff_meth.R

file.list <- list("recovery_0_1_CpG.txt","recovery_0_2_CpG.txt","recovery_0_4_CpG.txt",
                  "recovery_1_2_CpG.txt","recovery_2_1_CpG.txt","recovery_2_5_11_CpG.txt",
                  "recovery_2_5_9_CpG.txt", "recovery_3_5_15_CpG.txt","recovery_3_5_1_CpG.txt",
                  "recovery_3_5_2_CpG.txt", "recovery_3_6_CpG.txt",
                  "pesticide_6_2_CpG.txt","pesticide_6_3_CpG.txt","pesticide_7_3_CpG.txt",
                  "pesticide_7_5_CpG.txt","pesticide_7_5_4_CpG.txt","pesticide_8_5_3_CpG.txt",
                  "pesticide_9_5_1_CpG.txt", "pesticide_9_5_3_CpG.txt","pesticide_9_6_CpG.txt",
                  "pesticide_9_20_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("0_1", "0_2","0_4","1_2","2_1","2_5_11","2_5_9","3_5_15",
                                      "3_5_1","3_5_2","3_6",
                                      "6_2","6_3","7_3","7_5", "7_5_4","8_5_3","9_5_1",
                                      "9_5_3","9_6","9_20"),
                     treatment = c(rep(0,11), rep(1,10)),
                     assembly="daphmag_2.4",
                     context="CpG")

## -------------------------------------------------------------------------

# Filter data for outliers and coverage and normalise counts
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

normalized <- normalizeCoverage(filtered_data)

## -------------------------------------------------------------------------

# Only take CpGs present in all samples (maybe too stingent, can adjust later)
meth_all_data <- unite(normalized, destrand=TRUE) 
nrow(meth_all_data)
# Eu vs Pest = 236228
# Eu vs Pris = 144723
# Eu vs R = 108271
# Pris vs R = 101622
# R vs Pest = 162771
# Pest vs Pris = 203823

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
# This helps a LOT with the late FDR correction
df_meth_all <- getData(meth_all_data)

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
n <- df_meth_all[,41:42]
o <- df_meth_all[,44:45]
p <- df_meth_all[,47:48]
q <- df_meth_all[,50:51]
r <- df_meth_all[,53:54]
s <- df_meth_all[,56:57]
t <- df_meth_all[,59:60]
u <- df_meth_all[,62:63]
v <- df_meth_all[,65:66]# Check have right number of dataframes for number of sampels

bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,j,k,l,m,n,o,p,q,r,s,t,u,v)) { # Check have right number of dataframes for number of sampels
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions)
# Eu vs Pest = 5559 CpGs
# Eu vs Pris = 2431
# Eu vs R = 3871
# Pris vs R = 4338
# R vs Pest = 7331
# Pest vs Pris = 5661

subset_methBase <- methylKit::select(meth_all_data, meth_positions)

## -------------------------------------------------------------------------

# Save the dataframe for later use, including making nice plots
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="rec_vs_pest_objectmethbase.txt", quote=F, row.names = F, sep = '\t')

## -------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("rec_vs_pest_correlation.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("rec_vs_pest_cluster.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("rec_vs_pest_PCA.pdf")
PCASamples(subset_methBase)
dev.off()


## -------------------------------------------------------------------------

# Differential methylation (ensure to include any appropriate covariates here)
diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 1)
write.csv(diff_meth, file="rec_vs_pest_all_tested_meth_sites_MSCfilter.csv")

diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.05)
write.csv(diff_meth_10, file="rec_vs_pest_DMRs_min10percentDiff_qval0.05_MSCfilter.csv")
nrow(diff_meth_10)
# Eu vs Pest = 235 diff CpGs
# Eu vs Pris = 179
# Eu vs R = 162
# Pris vs R = 145
# R vs Pest = 239
# Pest vs Pris = 235


## -------------------------------------------------------------------------
# Annotating diff meth sites with gene IDs (can also include a bigger annotation file with exons etc. if you want more resolution)

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

