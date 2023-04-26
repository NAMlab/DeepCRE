# Deep learning the cis-regulatory code for gene expression in selected model plants
### Generating non-homologous validation set
These are the steps taken to generate non-homologous training and validation sets
training set.
1. Download protein fasta file from Ensembl Plants
```shell
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```
```shell
gunzip Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```
2. Create protein blast database and perform an all by all blast
```shell
makeblastdb -in Arabidopsis_thaliana.TAIR10.pep.all.fa -dbtype prot -title arabidopsis -parse_seqids -hash_index -out arabidopsis
```
```shell
blastp -db arabidopsis -query Arabidopsis_thaliana.TAIR10.pep.all.fa  -out Blast_ara_to_ara -outfmt 6
```
3. USE the python script `produce_non_homologous_val_sets.py` to create non-homologous validation set of gene_ids

### Mappping reads to reference transcriptomes using kallisto
One can use the helper script `map_reads.py` to automatically download fastq files from NCBI SRA and map to a reference
transcriptome. This script however assumes that user has created the index file with `kallisto index`.

Once mapping is down with kallisto, the helper r script `process_kallisto.R` can be used to produce a file with tpm 
counts across sra experiments

### Training CNN models and getting motifs
1. To train the CNN models for SSR and MSR please use `train_low_high_expressed_models.py` and 
`train_low_high_expressed_multi_species_models.py` respectively

2. Once models are trained, run `motif_discovery_ssr.py` and `motif_discovery_msr.py`to run deeplift
and generate motifs with modisco

3. run `extract_motifs_ssr.py` and `extract_motifs_msr.py` to get the motfs produced by modisco from
the h5 files, for further comparison with JASPAR motifs

### Training random forest classifiers
1. Build tabular dataset of generic features using `create_generic_features.py`
2. run `random_forest_ssr.py` and `random_forest_msr.py` to train binary and multi-taks classifiers