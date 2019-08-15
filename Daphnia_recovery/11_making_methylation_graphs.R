## -------------------------------------------------------------------------
## Making Fancy Genome-Methyation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetic_Recovery/Sequencing_Data")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Read in the subsetted file from methylkit (NEED TO COVER ALL 5 EXPERIMENTS)
objectmethbase1 <- read_delim("./Differential_meth/control_vs_control/C4_vs_C6_objectmethbase.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase1$chrBase <- paste(objectmethbase1$chr, ".", objectmethbase1$start, sep="")
objectmethbase1 <- objectmethbase1[,-3]
objectmethbase1$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase1$strand)


objectmethbase2 <- read_delim("./Differential_meth/pest_vs_pest/P4_vs_P6_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase2$chrBase <- paste(objectmethbase2$chr, ".", objectmethbase2$start, sep="")
objectmethbase2 <- objectmethbase2[,-3]
objectmethbase2$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase2$strand)


objectmethbase3 <- read_delim("./Differential_meth/control_vs_recovery/C6_vs_R6_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase3$chrBase <- paste(objectmethbase3$chr, ".", objectmethbase3$start, sep="")
objectmethbase3 <- objectmethbase3[,-3]
objectmethbase3$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase3$strand)


## -------------------------------------------------------------------------

# Setset out each sample (c4 and c6)
a <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
b <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
c <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
d <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
e <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
f <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
g <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs7","numTs7")]
h <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
j <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs9","numTs9")]
k <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
l <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
m <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]

# Setset out each sample (p4 and p6)
a1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
b1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
c1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
d1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
e1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
f1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
g1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs7","numTs7")]
h1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
j1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs9","numTs9")]
k1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
l1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
m1 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]

# Setset out each sample (r6, sample order P6,R6,P6,R6 ...)
a2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
b2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
c3 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
d2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
e2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
f2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]

## -------------------------------------------------------------------------

# Write out each sample
all_files <- list(a,b,c,d,e,f,g,h,j,k,l,m,a1,b1,c1,d1,e1,f1,g1,h1,j1,k1,l1,m1,a2,b2,c2,d2,e2,f2)

for(i in seq_along(all_files)){
  colnames(all_files[[i]])[c(3,5,6,7)] <- c("base","coverage","numCs","numTs")
  all_files[[i]]$freqC <- round((all_files[[i]]$numCs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]]$freqT <- round((all_files[[i]]$numTs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]] <- all_files[[i]][-c(6,7)]
  myfile <- file.path("./", paste0(i,"_","subsetted_final.txt"))
  write.table(all_files[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}

# Urgh need to go and rename each file with sample name and condition,
# check the methylkit scritps for the order of samples

## -------------------------------------------------------------------------
# Use methylkit to get the data all togther and plottable 
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetic_Recovery/Sequencing_Data/Genome_wide_meth/subsetted_inputs")

file.listA <- list("36_01_C4_subsetted_final.txt","36_01_C6_subsetted_final.txt",
                   "53_01_C4_subsetted_final.txt","53_01_C6_subsetted_final.txt",
                   "54_01_C4_subsetted_final.txt","54_01_C6_subsetted_final.txt",
                   "6_2_C4_subsetted_final.txt","6_2_C6_subsetted_final.txt",
                   "7_5_C4_subsetted_final.txt","7_5_C6_subsetted_final.txt",
                   "9_20_C4_subsetted_final.txt","9_20_C6_subsetted_final.txt",
                   "36_01_P4_subsetted_final.txt","36_01_P6_subsetted_final.txt",
                   "53_01_P4_subsetted_final.txt","53_01_P6_subsetted_final.txt",
                   "54_01_P4_subsetted_final.txt","54_01_P6_subsetted_final.txt",
                   "6_2_P4_subsetted_final.txt","6_2_P6_subsetted_final.txt",
                   "7_5_P4_subsetted_final.txt","7_5_P6_subsetted_final.txt",
                   "9_20_P4_subsetted_final.txt","9_20_P6_subsetted_final.txt",
                   "36_01_R6_subsetted_final.txt","53_01_R6_subsetted_final.txt",
                   "54_01_R6_subsetted_final.txt","6_2_R6_subsetted_final.txt",
                   "7_5_R6_subsetted_final.txt","9_20_R6_subsetted_final.txt")

sample_list <- list("36_01_C4","36_01_C6",
                    "53_01_C4","53_01_C6",
                    "54_01_C4","54_01_C6",
                    "6_2_C4","6_2_C6",
                    "7_5_C4","7_5_C6",
                    "9_20_C4","9_20_C6",
                    "36_01_P4","36_01_P6",
                    "53_01_P4","53_01_P6",
                    "54_01_P4","54_01_P6",
                    "6_2_P4","6_2_P6",
                    "7_5_P4","7_5_P6",
                    "9_20_P4","9_20_P6",
                    "36_01_R6","53_01_R6",
                    "54_01_R6","6_2_R6",
                    "7_5_R6","9_20_R6")

# Make a list of all genotypes in the right order = genotype_list

raw_data <- methRead(file.listA,
                     sample.id = sample_list,
                     treatment = c(rep(1:2, 6), rep(3:4, 6), rep(5, 6)),
                     assembly="bter_1.0", 
                     context="CpG")

meth_all_data <- unite(raw_data)#697 CpGs mmmm


## -------------------------------------------------------------------------
# Make a nice dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)
labs$Experiment <- labs$label
labs$Experiment <- gsub(".*_", "", labs$Experiment)

labs$Population<- c(rep("Pesticide", 7), "Pristine", rep("Pesticide", 5),
             "Pristine", rep("Pesticide", 2),rep("Pristine", 4),"Pesticide",
             rep("Pristine", 9))

ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=0.8,
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.1, colour=Experiment,fontface="bold"),
            show.legend = T, angle = 90, size = 6, hjust = 1)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))+
  ylim(-0.6,2.5)

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(meth_all_data, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)
PCA_data1$Population <- c(rep("Pristine", 6), rep("Pesticide",6),
                          rep("Pristine", 6), rep("Pesticide",6),
                          rep("Pristine", 3), rep("Pesticide",3))
PCA_data1$Experiment <- PCA_data1$sample
PCA_data1$Experiment <- gsub(".*_", "", PCA_data1$Experiment)
                          
                          
percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste( colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Population))+
  geom_point(size=5)+
  geom_text_repel(aes(label=sample), size=7,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25)+
  theme_bw()+
  xlab(paste0("PC1: ",percentage[1],"variance")) +
  ylab(paste0("PC2: ",percentage[2],"variance")) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank())


## -------------------------------------------------------------------------
# Look for overall differences in CpG methylation levels across experiments

meth <- read_delim("./Useful_Dmag_files/Dmag_methylated_genes_with_exp_and_pop_weighted_meth.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

# Remove rows where there are no values for a single experiment
meth_conservative <- meth[complete.cases(meth),] #3959 genes

library(reshape2)
meth_melted <- melt(meth_conservative)
colnames(meth_melted) <- c("geneID", "pop_exp","weighted_meth")

library(tidyr)
meth_metled1<- separate(data = meth_melted, col = pop_exp, 
                                   into = c("condition", "generation", "Population"), sep = "\\_")
meth_metled1$exp_gen <- paste(meth_metled1$condition, meth_metled1$generation, sep="_")

meth_metled1$Population[meth_metled1$Population == "pest"] <- "Pesticide"
meth_metled1$Population[meth_metled1$Population == "pris"] <- "Pristine"


ggplot(meth_metled1, aes(x=exp_gen, y=weighted_meth, fill=Population))+
  geom_boxplot()+
  xlab("Condition and Generation")+
  ylab("Mean Weighted Methylation Level per Gene")+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title=element_text(size = 18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+
  scale_x_discrete(breaks = c('control_4','control_6','pesticide_4', "pesticide_6", "recovery_6"),
                   labels = c('Control 4','Contol 6','Pesticide 4', "Pesticide 6", "Recovery 6"))

    

## HERE 
model1<-lm(weighted_meth ~ condition * generation * population, data=meth_metled1)
summary(model1) # no interactions

model2<-lm(weighted_meth ~ condition + generation + population, data=meth_metled1)
summary(model2)
# generation 6 is different ? :S

anova(model1, model2)


model3<-lm(weighted_meth ~ condition * population, data=meth_metled1)
summary(model3)

model4<-lm(weighted_meth ~ condition + population, data=meth_metled1)
summary(model4)

anova(model3, model4) # no interaction

summary.lm(model4) 
drop1(model4, test="F") # Recovery condition is different to control and pest :S

# ----------------------------------------------------------------------------------
# Make PCA and dendogram by weighted meth of genes not CpGs?







