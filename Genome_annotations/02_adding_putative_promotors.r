# ----------------------------------------------------------------
### Adding putative promotor information to the annotation file
# ----------------------------------------------------------------
setwd("~/Dropbox/Birmingham/DmagnaAnnotate26Nov2020")
library(readr)

annotation <- read_delim("Daphnia_magna_LRV0_1_longestIsoform_plusIntrons.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
tail(annotation)
annotation <- annotation[,-c(2,6,8)]
annotation <- annotation[!duplicated(annotation),] #0
colnames(annotation) <- c("chr","feature","start", "end","strand","gene_id")

upsteam5UTRs <- subset(annotation, feature == "five_prime_UTR")

# Double check for genes with two UTRs! Ahhh ok, it's when there are multiple TSS maybe?
look <- upsteam5UTRs[duplicated(upsteam5UTRs$gene_id),]
check <- annotation[annotation$gene_id %in% look$gene_id,]
range(table(check$gene_id[check$feature=="five_prime_UTR"])) # 2-15 UTRs in the duplicated ones (Nasonia) 2-9 for BB
check2 <- check[check$feature=="five_prime_UTR",]
length(unique(look$gene_id)) # 2484/12743

upsteam5UTRs_pos <- subset(upsteam5UTRs, strand == "+")
upsteam5UTRs_neg <- subset(upsteam5UTRs, strand == "-")

upsteam5UTRs_pos$promotor_start <- upsteam5UTRs_pos$start - 500
upsteam5UTRs_pos <- upsteam5UTRs_pos[upsteam5UTRs_pos$promotor_start > 0,]#0 genes (promotor overlaps scaffold start: removed)
upsteam5UTRs_pos <- upsteam5UTRs_pos[,-4]
colnames(upsteam5UTRs_pos)[3] <- "end"
colnames(upsteam5UTRs_pos)[6] <- "start"
upsteam5UTRs_pos$feature <- "promoter"

upsteam5UTRs_neg$promotor_start <- upsteam5UTRs_neg$end + 500
upsteam5UTRs_neg <- upsteam5UTRs_neg[upsteam5UTRs_neg$promotor_start > 0,]
upsteam5UTRs_neg <- upsteam5UTRs_neg[,-3] # rm redundent end column
colnames(upsteam5UTRs_neg)[3] <- "start"
colnames(upsteam5UTRs_neg)[6] <- "end"
upsteam5UTRs_neg$feature <- "promoter"

promoters <- rbind(upsteam5UTRs_pos,upsteam5UTRs_neg)

all <- rbind(annotation,promoters)

look <- all[all$gene_id== "Dmagna028519",] # neg gene
look <- all[all$gene_id== "Dmagna000013",] # pos gene


write.table(all, file="Daphnia_magna_LRV0_1_longestIsoform_plusIntrons_plusPromoters.txt ", 
            row.names = F, col.names = T, sep = '\t', quote = F)
