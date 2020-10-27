#----------------------------------------------
# Post-alignment ATAC-Seq QC
#----------------------------------------------

#library(BiocManager)
#BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments"))
#library(ATACseqQC)

setwd("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/ATAC/testing")


#----------------------------------------------
# Try something else
# FROM: https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html

library(magrittr)
library(dplyr)
library(ggplot2)
library(GenomicAlignments)

atacReads <- readGAlignmentPairs("0_1-1.bam", param = ScanBamParam(mapqFilter = 1, 
              flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
              "mapq", "isize")))
length(atacReads)
atacReads

atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)

fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
             Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
              geom_line() + theme_bw()
fragLenPlot

fragLenPlot + scale_y_continuous(trans = "log2")

fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
             437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
               theme_bw()

fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
              247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
             geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()


library(soGGi)

# Nucleosome free
nucFree <- regionPlot(bamFile = "0_1-1.bam", testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
                      forceFragment = 50)

# Mononucleosome
monoNuc <- regionPlot(bamFile = "0_1-1.bam", testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
                      forceFragment = 80)

# Dinucleosome
diNuc <- regionPlot(bamFile = "0_1-1.bam", testRanges = TSSs, style = "point", 
                    format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
                    forceFragment = 160)
