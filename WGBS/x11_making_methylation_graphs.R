## -------------------------------------------------------------------------
## Making Fancy Genome-Methyation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/Alignment_phased_oldref_Genome")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Read in the subsetted file from methylkit (NEED TO COVER ALL 4 POPS)
# Run this script for Pristine vs Eutrophic and for Pesticide vs Recovery to get all pops
objectmethbase1 <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/Alignment_phased_oldref_Genome/differential_meth_outputs/objectmethbase/eu_vs_pest_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase1$chrBase <- paste(objectmethbase1$chr, ".", objectmethbase1$start, sep="")
objectmethbase1 <- objectmethbase1[,-3]
objectmethbase1$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase1$strand)


objectmethbase2 <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/Alignment_phased_oldref_Genome/differential_meth_outputs/objectmethbase/eu_vs_pris_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase2$chrBase <- paste(objectmethbase2$chr, ".", objectmethbase2$start, sep="")
objectmethbase2 <- objectmethbase2[,-3]
objectmethbase2$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase2$strand)


objectmethbase3 <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/Alignment_phased_oldref_Genome/differential_meth_outputs/objectmethbase/eu_vs_rec_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase3$chrBase <- paste(objectmethbase3$chr, ".", objectmethbase3$start, sep="")
objectmethbase3 <- objectmethbase3[,-3]
objectmethbase3$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase3$strand)



## -------------------------------------------------------------------------

# Setset out each sample (eutrophic and pest)
sample_12_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
sample_12_4 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
sample_12_5_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
sample_13_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
sample_13_2 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
sample_13_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
sample_13_5_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs7","numTs7")]
sample_14_5_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
sample_15_5_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs9","numTs9")]

sample_6_2 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
sample_6_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
sample_7_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]
sample_7_5 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs13","numTs13")]
sample_7_5_4 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs14","numTs14")]
sample_8_5_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs15","numTs15")]
sample_9_5_1 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs16","numTs16")]
sample_9_5_3 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs17","numTs17")]
sample_9_6 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs18","numTs18")]
sample_9_20 <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs19","numTs19")]

# pristine
sample_36_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
sample_36_02 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
sample_48_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]
sample_48_02 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs13","numTs13")]
sample_53_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs14","numTs14")]
sample_54_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs15","numTs15")]
sample_54_02 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs16","numTs16")]
sample_74_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs17","numTs17")]
sample_77_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs18","numTs18")]
sample_88_01 <- objectmethbase2[,c("chrBase","chr","start","strand","coverage1","numCs19","numTs19")]

# recovery
sample_0_1 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
sample_0_2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
sample_0_4 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]
sample_1_2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs13","numTs13")]
sample_2_1 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs14","numTs14")]
sample_2_5_11 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs15","numTs15")]
sample_2_5_9 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs16","numTs16")]
sample_3_5_15 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs17","numTs17")]
sample_3_5_1 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs18","numTs18")]
sample_3_5_2 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs19","numTs19")]
sample_3_6 <- objectmethbase3[,c("chrBase","chr","start","strand","coverage1","numCs20","numTs20")]

## -------------------------------------------------------------------------

# Write out each sample
all_files <- mget(ls(pattern = "sample*"))

for(i in seq_along(all_files)){
  colnames(all_files[[i]])[c(3,5,6,7)] <- c("base","coverage","numCs","numTs")
  all_files[[i]]$freqC <- round((all_files[[i]]$numCs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]]$freqT <- round((all_files[[i]]$numTs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]] <- all_files[[i]][-c(6,7)]
  myfile <- file.path("./", paste0(names(all_files[i]),"_","subsetted_final.txt"))
  write.table(all_files[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}


## -------------------------------------------------------------------------
# Use methylkit to get the data all togther and plottable 

file.listA <- list("sample_12_3_subsetted_final.txt", "sample_12_4_subsetted_final.txt",	 "sample_12_5_1_subsetted_final.txt",  
"sample_13_1_subsetted_final.txt",	"sample_13_2_subsetted_final.txt",	"sample_13_3_subsetted_final.txt",	
"sample_13_5_1_subsetted_final.txt", "sample_14_5_1_subsetted_final.txt", "sample_15_5_1_subsetted_final.txt",
"sample_6_2_subsetted_final.txt","sample_6_3_subsetted_final.txt","sample_7_3_subsetted_final.txt",
"sample_7_5_subsetted_final.txt","sample_7_5_4_subsetted_final.txt","sample_8_5_3_subsetted_final.txt",
"sample_9_5_1_subsetted_final.txt", "sample_9_5_3_subsetted_final.txt","sample_9_6_subsetted_final.txt",
"sample_9_20_subsetted_final.txt",
"sample_36_01_subsetted_final.txt","sample_36_02_subsetted_final.txt","sample_48_01_subsetted_final.txt",
"sample_48_02_subsetted_final.txt","sample_53_01_subsetted_final.txt","sample_54_01_subsetted_final.txt",
"sample_54_02_subsetted_final.txt", "sample_74_01_subsetted_final.txt","sample_77_01_subsetted_final.txt",
"sample_88_01_subsetted_final.txt",
"sample_0_1_subsetted_final.txt","sample_0_2_subsetted_final.txt","sample_0_4_subsetted_final.txt",
"sample_1_2_subsetted_final.txt","sample_2_1_subsetted_final.txt","sample_2_5_11_subsetted_final.txt",
"sample_2_5_9_subsetted_final.txt", "sample_3_5_15_subsetted_final.txt","sample_3_5_1_subsetted_final.txt",
"sample_3_5_2_subsetted_final.txt", "sample_3_6_subsetted_final.txt")


raw_data <- methRead(file.listA,
                     sample.id = list("12_3", "12_4", "12_5_1", "13_1", "13_2", "13_3",
                                      "13_5_1", "14_5_1", "15_5_1",
                                      "6_2","6_3","7_3","7_5", "7_5_4","8_5_3","9_5_1",
                                      "9_5_3","9_6","9_20",
                                      "36_01","36_02","48_01","48_02", "53_01","54_01","54_02",
                                      "74_02","77_01","88_01",
                                      "0_1", "0_2","0_4","1_2","2_1","2_5_11","2_5_9","3_5_15",
                                      "3_5_1","3_5_2","3_6"),
                     treatment = c(rep(0, 9), rep(1, 10), rep(2, 10), rep(3, 11)),
                     assembly="daphmag_2.4",
                     context="CpG")

meth_all_data <- unite(raw_data)


## -------------------------------------------------------------------------
# Make a nice dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)
labs$Population <- c("Pesticide","Pristine","Pesticide","Pristine","Eutrophic","Eutrophic",
                     "Pesticide","Pristine","Eutrophic","Pristine","Pristine","Pristine","Pristine",
                     "Recovery","Recovery","Recovery","Eutrophic","Pristine","Pesticide","Pristine","Eutrophic",
                     "Pesticide","Recovery","Eutrophic","Eutrophic","Recovery","Recovery","Pesticide","Pesticide",
                     "Pesticide","Recovery","Recovery","Recovery","Pesticide","Recovery","Pesticide","Recovery",
                     "Eutrophic","Eutrophic","Pristine")

ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, size=0.8),
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.01, colour=Population, size=3, fontface="bold"),
            show.legend = F)+
  scale_color_manual(breaks = c("Pristine","Eutrophic","Pesticide","Recovery"),
                    values=c("#339966","#6600CC","#FFCC00", "#CC0066"))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.title = element_blank(),
        axis.text = element_blank())


## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(meth_all_data, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)
PCA_data1$Population <- c(rep("Eutrophic",9),rep("Pesticide",10),rep("Pristine",10),rep("Recovery",11))

percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste( colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Population))+
  geom_point(size=5)+
  geom_text_repel(aes(label=sample), size=6,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25)+
  scale_color_manual(breaks = c("Pristine","Eutrophic","Pesticide","Recovery"),
                     values=c("#339966","#6600CC","#FFCC00", "#CC0066"))+
  theme_bw()+
  xlab(paste0("PC1: ",percentage[1],"variance")) +
  ylab(paste0("PC2: ",percentage[2],"variance")) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank())


## -------------------------------------------------------------------------
# Look for overall differences in CpG methylation levels across pops

meth <- read_delim("Dmag_melted_weighted_meth_by_pop.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

# Remove rows where there are no values for a single pop
meth_conservative <- meth[complete.cases(meth),] #26010 genes

ggplot(meth_conservative, aes(x=origin, y=weightedMeth.mean, fill=origin))+
  geom_boxplot()+
  theme_bw()+
  xlab("Population")+
  ylab("Mean Weighted Methylation Level per Gene")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title=element_text(size = 18))+
  guides(fill=FALSE)+
  scale_fill_manual(values=c("#x","#x","#x", "#x"))


kruskal.test(plot_data, weightedMeth.mean ~ origin)

output= dunnTest( weightedMeth.mean ~ origin,
                  data=plot_data,
                  method="bh")
output


