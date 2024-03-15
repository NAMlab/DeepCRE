## Deep learning the cis-regulatory code for gene expression in selected model plants

[![DOI](https://zenodo.org/badge/632932657.svg)](https://zenodo.org/doi/10.5281/zenodo.10822013)

Please follow the steps below to reproduce the results from our work.
1. Download this repository.
2. Change directory into the **model** subdirectory, run the *fetch_genomes_and_annotation.sh* script. Firstly, this
will create 3 new subdirectories: genomes, gene_models, tpm_counts. Then, it will download genomes and gene models
for 4 plant species, uncompress them and store them in genomes and gene models subdirectories respectively.
3. Download expression counts for the project from supplementary data from publication and save these files 
within the tpm_counts subdirectory. NB: you should have 8 files, corresponding to 4 plant species and 2 tissues.
For example, for *Arabidopsis thaliana* you would have arabidopsis_counts.csv and arabidopsis_root_counts.csv 
for leaf and root tissues respectively.

### Training convolutional neural networks
- To train SSR and SSC models, run the *train_ssr_ssc_models_leaf.py* and *train_ssr_ssc_models_root.py* for
leaf and root tissues respectively.
- Train the MSR models using *train_msr_models_leaf.py* and *train_msr_models_root.py*.

Only after training CNN models can you run the scripts below that compute importance scores and generate motifs. Also 
note that deepLIFT which is used to compute importance scores is currently only compatible with tensorflow 1.x. So if
you build models with tensorflow 2.x, you won't be able to use these scripts.
### Computing importance scores and obtaining motifs
- Run the *motif_discovery_...* scripts for respective tissue and models.
- Then run *extract_motifs_ssr.py* or *extract_motifs_msr.py* to get the motifs out of the output produced by
modisco.

### Random forest models
- Firstly, create the features using the *create_generic_feature.py* script.
- Then run either *random_forest_msr.py* or *random_forest_ssr.py* for MSR and SSR models respectively.

### Investigate effect of different sequence lengths
To investigate the effects of different UTR or promoter sequence lengths, use the *effect_of_different_...*
scripts. These scripts will build several models based on different length specified within the scripts.


#### Generating validation_genes.pickle file
This file contains information of genes that have homologs only within their chromosomes, such that when we use
chromosome level cross validation for training, we mitigate the effects of homologs leaking information between our
training and test set. While I provide the pickle file used for this project, one can generate this themselves  by
firstly going into the data directory that sits as a sibling directory to model directory, then running the commands
below in the terminal:

```shell
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```
```shell
gunzip Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```
```shell
makeblastdb -in Arabidopsis_thaliana.TAIR10.pep.all.fa -dbtype prot -title arabidopsis -parse_seqids -hash_index -out arabidopsis
```
```shell
blastp -db arabidopsis -query Arabidopsis_thaliana.TAIR10.pep.all.fa  -out Blast_ara_to_ara -outfmt 6
```
The above assume that you have blast installed on your computer.
The above 4 lines are just for *Arabidopsis thaliana* but should be edited and repeated for the other 4 species. Once 
this is done, run the *produce_non_homologous_val_sets.py* script.
