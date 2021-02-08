# ----------------------------------------------------------------
### Adding putative promotor information to the annotation file
# ----------------------------------------------------------------
setwd("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020")
library(readr)

Daphnia_magna_LRV0_1 <- read_delim("Daphnia_magna_LRV0_1_longestIsoform_plusIntrons.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
head(Daphnia_magna_LRV0_1)
Daphnia_magna_LRV0_1 <- Daphnia_magna_LRV0_1[,-c(2,6,8)]
Daphnia_magna_LRV0_1 <- Daphnia_magna_LRV0_1[!duplicated(Daphnia_magna_LRV0_1),] #0
colnames(Daphnia_magna_LRV0_1) <- c("chr","feature","start", "end","strand","gene_id")

upsteam5UTRs <- subset(Daphnia_magna_LRV0_1, feature == "five_prime_UTR")
upsteam5UTRs <- subset(upsteam5UTRs, strand == "+")
upsteam3UTRs <- subset(Daphnia_magna_LRV0_1, feature == "three_prime_UTR")
upsteam3UTRs <- subset(upsteam3UTRs, strand == "-")

upsteam5UTRs$promotor_start <- upsteam5UTRs$start - 500
upsteam5UTRs <- upsteam5UTRs[upsteam5UTRs$promotor_start > 0,] #0 genes (promotor overlaps scaffold start: removed)

upsteam3UTRs$promotor_start <- upsteam3UTRs$start - 500
upsteam3UTRs <- upsteam3UTRs[upsteam3UTRs$promotor_start > 0,]

promoters <- rbind(upsteam3UTRs,upsteam5UTRs)

promoters <- promoters[,-4] # rm redundent end column
colnames(promoters)[3] <- "end"
colnames(promoters)[6] <- "start"
promoters$feature <- "promoter"

all <- rbind(Daphnia_magna_LRV0_1,promoters)

write.table(all, file="Daphnia_magna_LRV0_1_longestIsoform_plusIntrons_plusPromoters.txt ", 
            row.names = F, col.names = T, sep = '\t', quote = F)

# ----
# Add in the TE annotation as well
#TEs <- read_delim("Diaci_v3.0.ref.fa.mod.EDTA.intact.txt", 
#                  "\t", escape_double = FALSE, col_names = FALSE, 
#                  trim_ws = TRUE)

#head(TEs)
#TEs <- TEs[,-c(2,6,8)]
#TEs <- TEs[!duplicated(TEs),] #2
#TEs$X3 <- "TE"
#colnames(TEs) <- c("chr","feature","start", "end","strand","gene_id")

#all2 <- rbind(all, TEs)
#write.table(all2, file="Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters_plusTEs.txt ", 
#            row.names = F, col.names = T, sep = '\t', quote = F)
