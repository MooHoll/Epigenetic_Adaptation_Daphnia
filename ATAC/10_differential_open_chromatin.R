#------------------------------------------------------------------
# Differential open chromatin 
#------------------------------------------------------------------

setwd("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/ATAC/3_diff_chromatin")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(apeglm)
library(tidyr)
library(readr)

# Read in sample matrix and count data matrix
ATAC_sample_matrix <- read_delim("ATAC_sample_matrix.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_sample_matrix$tech_replicate[is.na(ATAC_sample_matrix$tech_replicate)] <- 1
ATAC_sample_matrix$collapse_by <- paste0(ATAC_sample_matrix$condition, "-", ATAC_sample_matrix$genotype, "-", ATAC_sample_matrix$bio_replicate)

ATAC_count_matrix <- read_delim("ATAC_count_matrix.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
ATAC_count_matrix <- as.data.frame(ATAC_count_matrix)
rownames(ATAC_count_matrix) <- ATAC_count_matrix$peak
ATAC_count_matrix <- ATAC_count_matrix[,-1]

# Convert 0 to 1 to avoid complete separation
ATAC_count_matrix[ATAC_count_matrix == 0] <- 1

# Make deseq2 object
dds= DESeqDataSetFromMatrix(countData = ATAC_count_matrix, colData = ATAC_sample_matrix, 
                            design = ~ population)
colData(dds) # 122 samples
colnames(dds)

# Collapse technical replicates
ddsColl <- collapseReplicates(dds, groupby= dds$collapse_by)
colData(ddsColl) # 79 samples

#------------------------------------------------------------------
# remove features with low counts - 79 samples, so say 
dds = ddsColl[ rowMeans(counts(ddsColl)) > 100, ] 
nrow(dds) # 58450/58568

# rlog transform counts (by average length and correcting for library size)
rld = rlog(dds, blind=FALSE)

#------------------------------------------------------------------
# PCA plot

# Modify the function to return PC3 and PC4 as well
# from: http://seqanswers.com/forums/showthread.php?t=66769
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


data = plotPCA(rld, intgroup = c("population"), ntop = 58450, returnData=TRUE)

percentVar = round(100 * attr(data, "percentVar"))
nrow(data)

# PC3 and PC4 also show a cloud
ggplot(data, aes(PC1, PC2, color=population)) + geom_point(size=10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed(ratio=2,clip = "on")+
  #geom_text_repel(aes(label=name), size=8,show.legend=FALSE, 
  #                point.padding = 2, box.padding = 1,
  #                segment.color = 'transparent') +
  scale_colour_manual("", breaks=c("CWP","EP","PP","PR"),
                      values = c("#DDCC77","#44AA99","#AA4499","#6699CC"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20))

# two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

# estimate size factors = normalize for dispersion
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$population, rld$sample, sep = " - " )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = F,
         fontsize = 5)

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$population, rld$sample, sep = " - " )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = F,
         fontsize = 5)

#------------------------------------------------------------------
# differential expression (used Benjamini-Hochberg adjustment)
dds = DESeq2::DESeq(dds, parallel=TRUE)

#------------------------------------------------------------------
# change the comparison (first group is more open if positive logFC)
res=results(dds, contrast=c("population", "PP", "PR"), alpha=0.05)
summary(res)

#distribution of coefficents of the model
plotMA(res, ylim=c(-4,4),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change', main="PP vs PR")


# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()

write.table(as.data.frame(res_ordered), file="diff_chrom_peaks_log2FC_PPvsPR.txt",
            sep="\t", quote = F, col.names = T, row.names = F)

#out of 58450 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #X total differentially open chromatin
nrow(res_significant[res_significant$log2FoldChange > 1.5,]) #X upregulated
nrow(res_significant[res_significant$log2FoldChange < -1.5,]) #X downregulated

#------------------------------------------------------------------
# Volcano plot
res_df <- as.data.frame(res)
res_df$padj[is.na(res_df$padj)]<-1
res_df$gene <- row.names(res_df)

res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)

res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-5,5)+
  ylim(0,6)+
  theme_bw()+
  ggtitle("PP vs PR")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none",
        title = element_text(size=20))


# Do some stats on the number of diff open chromatin regions
#Comparison Total	More open in first group	More open in second group
#CWP vs EP	257	183	74: X-squared = 46.23, df = 1, p-value = 1.052e-11
#CWP vs PP	0	0	0
#CWP vs PR	125	0	125: X-squared = 125, df = 1, p-value < 2.2e-16
#EP vs PP	0	0	0
#EP vs PR	107	0	107: X-squared = 107, df = 1, p-value < 2.2e-16
#PP vs PR	339	335	4:X-squared = 323.19, df = 1, p-value < 2.2e-16

# Goodness of fit
observed = c(335, 4)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected)
