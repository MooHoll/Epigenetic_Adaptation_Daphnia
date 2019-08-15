## -------------------------------------------------------------------------
## GO term enrichment of diff meth genes against all meth genes as background
## -------------------------------------------------------------------------

# NOTE: following Alun's script, give citation to Bebane et al. (2019)

# Load packages etc.
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/epigenetics/Epigenetics_Through_Time/GO_analysis")
library(GOstats)
library(GSEABase)
library(doBy)

## -------------------------------------------------------------------------

# Reading in background GO set 
GO_annotations <- read.table("Dmag_GOs_all_methylated_genes_as_background.txt")

GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]


## -------------------------------------------------------------------------

# Create a GO set object, gene set collection and universe

GO_frame <- GOFrame(GO_annotations,organism = "Daphnia magna")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))


## -------------------------------------------------------------------------

# Get genes for enrichment test

# Change file here depending on gene list
my_genes <- read.csv("x", header=T)
my_genes <- read.csv("x", header=T)
my_genes <- read.csv("x", header=T)

my_genes <- as.data.frame(na.omit(my_genes$geneID))
colnames(my_genes) <- "genes"

my_genes <- as.vector(my_genes[,1])


# Keep only annotated genes
my_genes <- my_genes[my_genes %in% universe]

# Pritine -> Eutrophic X/X
# Eutrophic -> Pesticide X/X
# Pesticide -> Recovery X/X


## -------------------------------------------------------------------------

# Parameters for hypergeometric test
Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Bumblebee Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      
      
      
      param_list <- c(param_list,parameters)
      
      
    }
  }
  
  names(param_list) <- name_1
  return(param_list)
}
param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


## -------------------------------------------------------------------------

# Carry out hypergeometric test

Hyper_G_test <- function(param_list){
  
  Hyper_G_list <- list()
  
  for(i in 1:length(param_list)){
    
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
    
  }
  
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


## -------------------------------------------------------------------------

# FDR correction and get enriched biological process terms
Result <- summary(GO_enrichment[["BP_over"]])

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

REVIGO <- Result_FDR[,1:2]

write.table(REVIGO,"pristine_eutrophic_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"eutrophic_pesticide_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"pesticide_recovery_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)


# NOTE: can use REVIGO or go back and merge with the original GO file which contains the annotation of the GO ids
