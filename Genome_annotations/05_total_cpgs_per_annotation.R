## -------------------------------------------------------------------------
## Weighted Methylation per annotation: make file with cpg count per annoation
## -------------------------------------------------------------------------

# Load packages etc.
library(sqldf)
library(readr)
library(doBy)

## -------------------------------------------------------------------------
# Making the file which has the total CpGs per gene information

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("chr", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)
head(cpgs)

# --------------------------------------------------------------------

genes <- read.delim("Daphnia_magna_LRV0_1_ALL_ANNOTATIONS.txt", header=T)
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)
head(genes)

output <- sqldf("SELECT cpgs.chr,
                cpgs.cpg_position,
                genes.chr,
                genes.feature,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.chr = genes.chr
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ gene_id+chr+start+end+feature, data = output, FUN=sum)

write.table(final, file="Daphnia_magna_LRV0_1_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)
