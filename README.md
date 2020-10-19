# The potential for epigenetic adaptation in the waterflea, *Daphnia magna*.

This project is part of a large collaborative effort run by Dr. Luisa Orsini and Professor John Colbourne at the University of Birmingham, UK. This repository holds all code and analyses used to attempt to identify epigenetic mechanisms of adaptaion. 

---

## NOTES:
We now have a new version of the *D. magna* genome at chromosome level (10 chromosomes + 61 scaffolds). I will therefore re-run all current analyses using this version. In order to progress further we need the annotation, which Vignesh is currently working on. 

---

## The Plan Methylation Data:
- Genome-wide DNA methylation differences between populations (all combinations)
    - Overall levels including non-CpG (boxplots)
    - Chromosome level differences (including sex chromosomes) (I like Mather's et al. density plot better than our current one)
    - PCA clustering
    - Mean feature levels (e.g. promotors/exons etc.) **Before doing this we want to add in TE annotation probably using this new pipeline [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y)**
    - Counts in high/medium/low/no bins per feature
    - Exon/intron average levels per gene **code needs writing from scratch**
- Differential DNA methylation on the CpG and feature level
    - Location of differentially methylated CpGs and count
    - Relative position in gene (if enriched there)
    - MA / scatter plot 
    - Gene ontology (GO) enrichment of differentially methylated genes **Need new GO terms from new annotation file when ready**

---

## The Plan ATAC Data:
- Differential chromatin profiles between populations
    - Find paper that's done this well and check out graphs etc.

---

## The Plan Integration Meth and Genome Data:
- 2/3rd codon useage and prevalence of DNA methylation if CpG
- More generally is DNA meth associated with SNPs (we have the SNPs but how to figure this out...)

---

## The Plan Integration ATAC and Genome Data:
- This analysis will be a version of above

---

## The Plan Integration Meth, ATAC and Genome Data:
- Oh my days no idea

---

NOTE: we are currently not planning to add gene expression or phenotypic data (such as mobility), but we may do this depending on the results.

---

For more information please contact:<br/>
hollie_marshall@hotmail.co.uk
