#--------------------------------------------------------------------
# Making intergenic regions or unannotated regions to quantify background methylation
#--------------------------------------------------------------------

setwd("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020")
library(readr)
library(dplyr)

merged_annotations <- read_delim("Daphnia_magna_LRV0_1_longestIsoform_plusIntrons_plusPromoters_plusTEs.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
head(merged_annotations)

# We can remove some annotations which are within others to decrease the file size
annotations_cleaner <- merged_annotations[merged_annotations$feature == "gene"|
                                            merged_annotations$feature =="promoter",]
annotations_cleaner <- annotations_cleaner[,-c(2,5)]
              
# Order the file and make in a format the python script wants, from promoter to 3' UTR
ordered <- annotations_cleaner %>% arrange(chr, start)
head(ordered)
colnames(ordered) <- c("chr","start","end","gene")
ordered <- ordered[,c(4,1,2,3)]

write.table(ordered, file="annotations_for_getting_intergenic_regions.gtf",
            col.names = T, row.names = F, quote = F, sep = '\t')

# Run getting_intergenic_regions.py as:
# python 03_getting_intergenic_regions.py annotations_for_getting_intergenic_regions.gtf > intergenic.bed

intergenic <- read_delim("intergenic.bed", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

# Remove rows which are incorrect becasue they wer created when a TE is found
# within gene conordinates
head(intergenic)

intergenic$remove <- ifelse(intergenic$X2 > intergenic$X3, "yes","no")
intergenic$remove_aswell <- ifelse(intergenic$X2 == intergenic$X3, "yes","no")
intergenic_for_sure <- intergenic[intergenic$remove == "no",]
intergenic_for_sure <- intergenic_for_sure[intergenic_for_sure$remove_aswell == "no",]

# Neaten the file up for later use
intergenic_for_sure$X4 <- "intergenic"
intergenic_for_sure <- intergenic_for_sure[,-c(5,6)]
colnames(intergenic_for_sure) <- c("chr","start","end","feature")
intergenic_for_sure$gene_id <- "intergenic"

# add to other annotations
head(merged_annotations)
intergenic_for_sure$strand <- NA

all <- rbind(merged_annotations, intergenic_for_sure)

write.table(all, file="Daphnia_magna_LRV0_1_ALL_ANNOTATIONS.txt",
            col.names = T, row.names = F, quote = F, sep = '\t')
