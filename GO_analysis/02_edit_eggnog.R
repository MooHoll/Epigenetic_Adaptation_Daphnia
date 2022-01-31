# Make a useable file from the excel of the eggNOG output

# NOTE: First removed all info from the .xlxs file except gene id and GO annotation

setwd("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/GO_analysis")

library(tidyr)
library(dplyr)
library(readxl)

out_emapper_annotations_copy <- read_excel("~/Dropbox/Birmingham/1_Epigenetics_Through_Time/2020_analyses/GO_analysis/out.emapper.annotations.xlsx", 
                                           col_names = TRUE)
colnames(out_emapper_annotations_copy) <- c("gene_id","GO_term")

new <- out_emapper_annotations_copy %>% 
  mutate(GO_term = strsplit(as.character(GO_term), ",")) %>% 
  unnest(GO_term)

new <- new[!new$GO_term=="-",]

write.table(new, file="dmagna_GO_annotations.txt", sep="\t", quote = F, row.names = F,
            col.names = T)
