# -------------------------------------------------------------------------
# Using Christian's version of this script to pull out methylated bases for pop epigenetics
# -------------------------------------------------------------------------

# On alice: qsub -I -l walltime=01:00:00 -l pvmem=100gb -l nodes=1:ppn=1

#Load Packages
library(methylKit)
library(readr)
library(dplyr)

# -------------------------------------------------------------------------

file.list <- list("eutrophic-12_3_CpG.txt","eutrophic-12_4_CpG.txt","eutrophic-12.5_1_CpG.txt",
                  "eutrophic-13_1_CpG.txt", "eutrophic-13_2_CpG.txt", "eutrophic-13_3_CpG.txt",
                  "eutrophic-13.5_1_CpG.txt", "eutrophic-14.5_1_CpG.txt", "eutrophic-15.5_1_CpG.txt",
                  "pesticide-6_2_CpG.txt", "pesticide-6_3_CpG.txt", "pesticide-7_3_CpG.txt",
                  "pesticide-7.5_4_CpG.txt", "pesticide-7_5_CpG.txt", "pesticide-8.5_3_CpG.txt",
                  "pesticide-9_20_CpG.txt", "pesticide-9.5_1_CpG.txt", "pesticide-9.5_3_CpG.txt",
                  "pesticide-9_6_CpG.txt", "pristine-36_01_CpG.txt", "pristine-36_02_CpG.txt",
                  "pristine-48_01_CpG.txt", "pristine-48_02_CpG.txt", "pristine-53_01_CpG.txt",
                  "pristine-54_01_CpG.txt", "pristine-54_02_CpG.txt", "pristine-74_01_CpG.txt",
                  "pristine-77_01_CpG.txt", "pristine-88_01_CpG.txt", "recovery-0_1_CpG.txt",
                  "recovery-0_2_CpG.txt", "recovery-0_4_CpG.txt", "recovery-1_2_CpG.txt",
                  "recovery-2_1_CpG.txt", "recovery-2.5_11_CpG.txt", "recovery-2.5_9_CpG.txt",
                  "recovery-3.5_15_CpG.txt", "recovery-3.5_1_CpG.txt", "recovery-3.5_2_CpG.txt",
                  "recovery-3_6_CpG.txt")

treat_conditions <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3)

raw_data <- methRead(file.list,
                     sample.id = list("EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                                      "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                                      "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                                      "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                                      "PP_LRV7_5_4", "PP_LRV7_4", "PP_LRV8_5_3",
                                      "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                                      "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                                      "PR_LR2_48_01", "PR_LR2_48_02", "PR_LR3_53_01",
                                      "PR_LR2_54_01", "PR_LR2_54_02", "PR_LR3_74_01",
                                      "PR_LR3_77_01", "PR_LR3_88_01", "CWP_LRV0_1",
                                      "CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                                      "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                                      "CWP_LRV3_5_15", "CWP_LRV3_5_1", "CWP_LRV3_5_2",
                                      "CWP_LRV3_6"),
                     treatment = treat_conditions,
                     assembly="daphmag_0_1",
                     context="CpG")

filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
rm(raw_data)

normalized <- normalizeCoverage(filtered_data)
rm(filtered_data)

# Remove min per group so we keep all data
meth_all_data <- unite(normalized, destrand=TRUE, min.per.group = 0L)

df_meth_all <- getData(meth_all_data) 
nrow(df_meth_all) # for all 4880224 positions identified in the genome as CpG positions
#write.table(df_meth_all, file="all_raw_data.txt", sep="\t", quote=F, col.names = T, row.names = T)
rm(meth_all_data)

e1 <- df_meth_all[,c(1:2,5:6)]
e2 <- df_meth_all[,c(1:2,8:9)]
e3 <- df_meth_all[,c(1:2,11:12)]
e4 <- df_meth_all[,c(1:2,14:15)]
e5 <- df_meth_all[,c(1:2,17:18)]
e6 <- df_meth_all[,c(1:2,20:21)]
e7 <- df_meth_all[,c(1:2,23:24)]
e8 <- df_meth_all[,c(1:2,26:27)]
e9 <- df_meth_all[,c(1:2,29:30)]
pe1 <- df_meth_all[,c(1:2,32:33)]
pe2 <- df_meth_all[,c(1:2,35:36)]
pe3 <- df_meth_all[,c(1:2,38:39)]
pe4 <- df_meth_all[,c(1:2,41:42)]
pe5 <- df_meth_all[,c(1:2,44:45)]
pe6 <- df_meth_all[,c(1:2,47:48)]
pe7 <- df_meth_all[,c(1:2,50:51)]
pe8 <- df_meth_all[,c(1:2,53:54)]
pe9 <- df_meth_all[,c(1:2,56:57)]
pe10 <- df_meth_all[,c(1:2,59:60)]
pr1 <- df_meth_all[,c(1:2,62:63)]
pr2 <- df_meth_all[,c(1:2,65:66)]
pr3 <- df_meth_all[,c(1:2,68:69)]
pr4 <- df_meth_all[,c(1:2,71:72)]
pr5 <- df_meth_all[,c(1:2,74:75)]
pr6 <- df_meth_all[,c(1:2,77:78)]
pr7 <- df_meth_all[,c(1:2,80:81)]
pr8 <- df_meth_all[,c(1:2,83:84)]
pr9 <- df_meth_all[,c(1:2,86:87)]
pr10 <- df_meth_all[,c(1:2,89:90)]
r1 <- df_meth_all[,c(1:2,92:93)]
r2 <- df_meth_all[,c(1:2,95:96)]
r3 <- df_meth_all[,c(1:2,98:99)]
r4 <- df_meth_all[,c(1:2,101:102)]
r5 <- df_meth_all[,c(1:2,104:105)]
r6 <- df_meth_all[,c(1:2,107:108)]
r7 <- df_meth_all[,c(1:2,110:111)]
r8 <- df_meth_all[,c(1:2,113:114)]
r9 <- df_meth_all[,c(1:2,116:117)]
r10 <- df_meth_all[,c(1:2,119:120)]
r11 <- df_meth_all[,c(1:2,122:123)]

datalist = list(e1,e2,e3,e4,e5,e6,e7,e8,e9,pe1,pe2,pe3,pe4,pe5,pe6,pe7,pe8,pe9,pe10,pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11)

sample.id = list("EP_LRV12_3","EP_LRV12_4","EP_LRV12_5_1", 
                 "EP_LRV13_1", "EP_LRV13_2", "EP_LRV13_3",
                 "EP_LRV13_5_1", "EP_LRV14_5_1", "EP_LRV15_5_1",
                 "PP_LRV6_2", "PP_LRV6_3", "PP_LRV7_3",
                 "PP_LRV7_5_4", "PP_LRV7_4", "PP_LRV8_5_3",
                 "PP_LRV9_20", "PP_LRV9_5_1", "PP_LRV9_5_3",
                 "PP_LRV9_6", "PR_LR2_36_01", "PR_LR2_36_02",
                 "PR_LR2_48_01", "PR_LR2_48_02", "PR_LR3_53_01",
                 "PR_LR2_54_01", "PR_LR2_54_02", "PR_LR3_74_01",
                 "PR_LR3_77_01", "PR_LR3_88_01", "CWP_LRV0_1",
                 "CWP_LRV0_2", "CWP_LRV0_4", "CWP_LRV1_2",
                 "CWP_LRV2_1", "CWP_LRV2_5_11", "CWP_LRV2_5_9",
                 "CWP_LRV3_5_15", "CWP_LRV3_5_1", "CWP_LRV3_5_2",
                 "CWP_LRV3_6")
names(datalist) <- sample.id

# Probability of success = average non-conversion rate
bt <- function(a, b, p = 0.0036) {binom.test(a, b, 0.0036, alternative="greater") $p.value}

for (i in seq_along(datalist)) {  
  colnames(datalist[[i]]) <- c("chr","position","coverage", "Ccount")
  datalist[[i]] <- datalist[[i]][!is.na(datalist[[i]]$coverage),]
  datalist[[i]]$pVal <- mapply(bt, datalist[[i]]$Ccount, datalist[[i]]$coverage)
  datalist[[i]]$FDR <- p.adjust(datalist[[i]]$pVal, method = "BH", n = length(datalist[[i]]$pVal))
  datalist[[i]]$methylated <- 0
  datalist[[i]]$methylated[datalist[[i]]$FDR < 0.05] <- 1
  datalist[[i]] <- datalist[[i]][,c(1,2,7)]
  colnames(datalist[[i]]) <- c("chr","position",paste0(names(datalist[i])))
}

all = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
write.table(all, file="methylation_calls_per_sample.txt", quote=F, sep="\t", row.names=F, col.names=T)

