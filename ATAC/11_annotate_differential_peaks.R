#------------------------------------------------------------------
# Label diff open chromatin regions with annotations
#------------------------------------------------------------------

setwd("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/ATAC/3_diff_chromatin")

library(readr)
library(sqldf)
library(dplyr)
library(tidyr)
library(ggplot2)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

genome_annotation<-read.csv.sql("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020/Daphnia_magna_LRV0_1_ALL_ANNOTATIONS.txt",
                                sql ="select * from file", sep="\t",header = T)
genome_annotation <- genome_annotation[!genome_annotation$feature=="CDS",]
genome_annotation <- genome_annotation[!genome_annotation$feature=="mRNA",]
genome_annotation <- genome_annotation[!genome_annotation$feature=="tRNA",]

files <- list.files(".", pattern="diff_chrom_peaks*")
samples <- lapply(files, function(x)read.table(x, sep="\t",header = T))
names(samples) <- c("CWPvsEP","CWPvsPP","CWPvsPR","EPvsPP","EPvsPR","PPvsPR")

for (i in seq_along(samples)) { 
  samples[[i]]$comparison <- paste0(names(samples[i]))
}

all <- bind_rows(samples)
significant <- subset(all, padj < 0.05)
significant <- significant[significant$log2FoldChange > 1.5 | significant$log2FoldChange < -1.5,]

significant <- separate(significant, gene, c("chr", "start", "end"),
                        sep = "-", fill = 'right')
significant <- significant[,c(2,7:10)]

# Check if any of the peaks overlap with annotated features
output_start <- sqldf("SELECT significant.chr,
                significant.start,
                significant.end,
                significant.comparison,
                significant.log2FoldChange,
                ga.chr,
                ga.start,
                ga.end,
                ga.feature,
                ga.gene_id
                FROM significant AS significant
                LEFT JOIN genome_annotation AS ga 
                ON significant.chr = ga.chr
                AND (significant.start >= ga.start AND significant.start <= ga.end)") 

output_end <- sqldf("SELECT significant.chr,
                significant.start,
                significant.end,
                significant.comparison,
                significant.log2FoldChange,
                ga.chr,
                ga.start,
                ga.end,
                ga.feature,
                ga.gene_id
                FROM significant AS significant
                LEFT JOIN genome_annotation AS ga 
                ON significant.chr = ga.chr
                AND (significant.start >= ga.start AND significant.start <= ga.end)")


head(output_start)
head(output_end)
output_start <- output_start[!is.na(output_start$gene_id),]
output_end <- output_end[!is.na(output_end$gene_id),]

both <- rbind(output_start, output_end)
both <- both[!duplicated(both),] 

# Just look how many out of 828 are annotated
annotated <- both[,c(1:4)]
annotated <- annotated[!duplicated(annotated),] # 726

head(both)
both <- both[,-c(1)]

both_plot <- both[!both$feature=="gene",] 
ggplot(both_plot, aes(x=feature, fill=comparison))+
  geom_bar(stat="count")+
  xlab("Genomic Feature")+
  ylab("Number of Significant Differential Chromatin Peaks")+
  theme_bw()+
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


write.table(both_plot, file="differntial_chromatin_peaks_annotated_with_features.txt",
            sep="\t", quote = F, col.names = T, row.names = F)

# Where are the differential peaks?
# What genes are they associated with?
# Are any features shared?


# Take gene and intergenic and plot on the chromosomes per comparison
for_plot <- both[(both$feature =="gene" | both$feature=="intergenic" | both$feature=="promoter"),]
head(for_plot) #728

for_plot_genes <- for_plot[for_plot$feature=="gene",]
for_plot_intergenic <- for_plot[for_plot$feature=="intergenic",]
#for_plot <- separate(for_plot, chr, c("scaffold", "number"),
#                        sep = "_", fill = 'right')
#for_plot_main_chrs <- for_plot[for_plot$number <=10,]
#for_plot_main_chrs$chr <- paste(for_plot_main_chrs$scaffold, for_plot_main_chrs$number, sep="_")


ggplot(for_plot, aes(fill=comparison))+
         geom_bar(aes(x = forcats::fct_infreq(chr)))+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Significant Differential Chromatin Peaks")+
  #ggtitle("Differential Chromatin Within Genes")+
  theme_bw()+
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())
  

