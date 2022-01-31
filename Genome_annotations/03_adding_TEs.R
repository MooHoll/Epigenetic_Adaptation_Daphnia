# ---------------------------------------------------------
# Adding TEs to the Daphnia full annotation file
# ---------------------------------------------------------
setwd("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020")
# From Vignesh's EDTA TE annotation
# cut -f 1,3,4,5,8,9 Daphnia_magna_LRV0_1.scaffolds.fa.mod.EDTA.intact.gff3 > te.txt
# sed 's/ID=.*Classification=//g' te.txt > te1.txt
# sed 's/;.*$//g' te1.txt > te2.txt

library(readr)

tes <- read_delim("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020/te2.txt", 
                  "\t", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE)

tes$X7 <- paste(tes$X2, tes$X6, sep=";")
tes$X2 <- "TE"
tes <- tes[,-6]
colnames(tes)<- c("chr","feature","start","end","strand","gene_id")

annotation <- read_delim("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020/Daphnia_magna_LRV0_1_longestIsoform_plusIntrons_plusPromoters.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

annotation_both <- rbind(annotation, tes)
write.table (annotation_both, file="Daphnia_magna_LRV0_1_longestIsoform_plusIntrons_plusPromoters_plusTEs.txt",
             sep="\t", quote = F, col.names = T, row.names = F)
