## -------------------------------------------------------------------------
## Making Fancy Genome-Methyation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("./")
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
objectmethbase1 <- read_delim("X_vs_Y_objectmethbase.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase$chrBase <- paste(objectmethbase$chr, ".", objectmethbase$start, sep="")
objectmethbase <- objectmethbase[,-3]
objectmethbase$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase$strand)


## -------------------------------------------------------------------------

# Setset out each sample
a <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
b <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
c <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
d <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
e <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
f <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
g <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs7","numTs7")]
h <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
i <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs9","numTs9")]
j <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs10","numTs10")]
k <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs11","numTs11")]
l <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs12","numTs12")]
m <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs13","numTs13")]
n <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs14","numTs14")]
o <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs15","numTs15")]
p <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs16","numTs16")]
q <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs17","numTs17")]
r <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs18","numTs18")]
s <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs19","numTs19")]
t <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs20","numTs20")]


## -------------------------------------------------------------------------

# Write out each sample
all_files <- list(letters[1-20])

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


## -------------------------------------------------------------------------
# Use methylkit to get the data all togther and plottable 

file.listA <- list("*final.txt")

# Make a list of all genotypes in the right order = genotype_list

raw_data <- methRead(file.listA,
                     sample.id = genotype_list,
                     treatment = c(rep(0, 10), rep(1, 10), rep(2, 10), rep(3, 10)),
                     assembly="bter_1.0", 
                     context="CpG")

meth_all_data <- unite(raw_data)


## -------------------------------------------------------------------------
# Make a nice dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)
labs$group <- c(rep("Recovery", 10), rep("Pesticide",10),rep("Eutrophic",10),rep("Pristine", 10))

labs$label<- genotype_list

ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, size=0.8),
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.01, colour=group, size=3, fontface="bold"),
            show.legend = F)+
  scale_color_manual(values=c("#x","#x","#x", "#x"))+
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
PCA_data1$Population <- c(rep("Recovery", 10), rep("Pesticide",10),rep("Eutrophic",10),rep("Pristine", 10))

percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste( colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Population))+
  geom_point(size=5)+
  geom_text_repel(aes(label=sample), size=6,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25)+
  scale_color_manual(values=c("#x","#x","#x", "#x"))+
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
meth_conservative <- meth_melted[complete.cases(meth_melted),] #X genes

ggplot(meth_conservative, aes(x=origin, y=weightedMeth.mean, fill=origin))+
  geom_boxplot()+
  theme_bw()+
  xlab("Population")+
  ylab("Mean Weighted Methylation Level per Gene")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title=element_text(size = 18))+
#  scale_x_discrete(labels = c('Male','Queen','Reproductive \n Worker'))+
  guides(fill=FALSE)+
  scale_fill_manual(values=c("#x","#x","#x", "#x"))


kruskal.test(plot_data, weightedMeth.mean ~ origin)

output= dunnTest( weightedMeth.mean ~ origin,
                  data=plot_data,
                  method="bh")
output


