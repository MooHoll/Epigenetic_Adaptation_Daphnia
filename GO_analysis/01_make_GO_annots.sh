# Make GO annotations

# on alice
conda create -n agat python=3.8
conda activate agat
conda install -c bioconda agat

# interactive job
agat_sp_extract_sequences.pl --gff Daphnia_magna_LRV0_1.gff3 -f Daphnia_magna_LRV0_1.scaffolds.fa -p -o dmagna_protein_seqs.fa

# Put protein seqs into eggnog mapper with standard options
# http://eggnog-mapper.embl.de

# This is the file we then want to download: out.emapper.annotations copy.xlsx