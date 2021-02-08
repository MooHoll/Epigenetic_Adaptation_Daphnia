# Making a bed-type file of the annotation for manipulation

# keep only the longest transcript (Edinburgh qm machine)
conda activate hmarshall 
conda install -c bioconda agat
agat_sp_keep_longest_isoform.pl -gff Daphnia_magna_LRV0_1.gff3 -o Daphnia_magna_LRV0_1_longestIsoformOnly.gff3
agat_sp_add_introns.pl -gff Daphnia_magna_LRV0_1_longestIsoformOnly.gff3 -o Daphnia_magna_LRV0_1_longestIsoform_plusIntrons.gff3

# make file to get intergenic regions
grep -v "#" Daphnia_magna_LRV0_1_longestIsoform_plusIntrons.gff3 > annotation1
sed 's/ID=.*Parent=//g' annotation1 > annotation2
sed 's/ID=//g' annotation2 > annotation3
sed 's/;.*$//g' annotation3 > annotation4
sed 's/-T.*$//g' annotation4 > Daphnia_magna_LRV0_1_longestIsoform_plusIntrons.txt

# VIGNESH IS MAKING THE TE ANNOTATION FILE CURRENTLY
# Fix the TE annotation file as well
#sed 's/;Sequence_ontology.*$//g' Diaci_v3.0.ref.fa.mod.EDTA.intact.gff3 > te_annot1
#sed 's/ID=.*Classification=//g' te_annot1 > Diaci_v3.0.ref.fa.mod.EDTA.intact.txt


