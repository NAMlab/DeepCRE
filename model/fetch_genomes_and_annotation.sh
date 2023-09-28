#! /bin/bash

# creating the required directories
echo "Creating required directories"
mkdir genomes
mkdir gene_models
mkdir tpm_counts

# Downloading reference genomes from Ensembl plants
echo "Downloading reference genomes"
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# Downloading corresponding annotation files in gtf format from Ensembl plants
echo "Downloading reference annotations"
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gtf/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.52.gtf.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gtf/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.52.gtf.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-52/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf.gz

gunzip *.gz
mv *.gtf gene_models
mv *toplevel.fa genomes

