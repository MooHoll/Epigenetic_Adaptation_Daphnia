## -------------------------------------------------------------------------
## Getting background methylated genes GO set
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/GO_analysis")
library(readr)

## -------------------------------------------------------------------------

# Read in overall GO file and methylated gene list
GO_annotations <- read_delim("Dmagna_allGOterms.txt", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
GO_annotations <- GO_annotations[,1:2]
colnames(GO_annotations) <- c("geneID","goID")

methylated_genes <- read_csv("Dmag_methylated_genes.txt")
colnames(methylated_genes) <- "geneID"


## -------------------------------------------------------------------------

# Merge to get GO terms associated with methylated genes

meth_gos <- merge(GO_annotations, methylated_genes, by="geneID")
length(meth_gos[is.na(meth_gos$geneID)]) 

write.table(meth_gos, file = "Dmag_GOs_all_methylated_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)


