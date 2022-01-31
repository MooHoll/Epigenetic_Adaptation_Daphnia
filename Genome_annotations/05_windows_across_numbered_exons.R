# ------------------------------------------------------
# Making that fancy exon and intron methylation graph
# ------------------------------------------------------

setwd("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020")

library(readr)
library(dplyr)
library(tidyr)

# ------------------------------------------------------
# check the filke is sorted so the 1st exon is 1st etc. 
# also check for strandedness as this could mess up the order
annotation <- read_delim("Daphnia_magna_LRV0_1_ALL_ANNOTATIONS.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
head(annotation)

# ------------------------------------------------------
# Take only the exons
exons_only <- subset(annotation, feature=="exon")

# Account for strandedness
positive_exons <- exons_only[exons_only$strand=='+',]
negative_exons <- exons_only[exons_only$strand=='-',]

# Number each exon per gene (positive file is already in order)
positve_numbered_exons <- positive_exons %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

# Re- order -ve strand exons so the last one goes first in the file order
negative_exons <- negative_exons %>% arrange(chr, -start)
head(negative_exons)
negative_numbered_exons <- negative_exons %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

look <- subset(negative_numbered_exons, gene_id =="Dmagna003752")

numbered_exons <- rbind(positve_numbered_exons, negative_numbered_exons)

# Label any exons >=6 as just n
numbered_exons$number[numbered_exons$number>=6] <- "n"

# ------------------------------------------------------
# Do the same for introns
introns_only <- subset(annotation, feature=="intron")

positive_introns <- introns_only[introns_only$strand=='+',]
negative_introns <- introns_only[introns_only$strand=='-',]

positve_numbered_introns <- positive_introns %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

negative_introns <- negative_introns %>% arrange(chr, -start)
negative_numbered_introns <- negative_introns %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

look <- subset(negative_numbered_introns, gene_id =="Dmagna003752")

numbered_introns <- rbind(positve_numbered_introns, negative_numbered_introns)

numbered_introns$number[numbered_introns$number>=6] <- "n"

# ------------------------------------------------------
# Replace the exons and introns in the main file

# Remove the original exons and introns
annotation1 <- annotation[!(annotation$feature == "exon" | annotation$feature=="intron"),]

# Put in the new ones
annotation1$number <- NA
final_annotation <- rbind(annotation1, numbered_exons, numbered_introns)

#Write out the file
write.table(final_annotation, file ="Daphnia_magna_LRV0_1_ALL_ANNOTATIONS_numbered_exons.txt", sep="\t",
            quote = F, col.names = T, row.names = F)

