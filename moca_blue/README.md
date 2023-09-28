# moca_blue
[2023-06-07]

MOCA BLUE

MOtif
  Characterization
&  Annotation
   from DEEP LEARNING feature enrichment

Welcome the the moca_blue suite!
from Simon M. Zumkeller
RStudio
2022.07.2 Build 576

This is a tool-box for the analyses of DNA motifs
that have been derived from deep-learning model features extraction.
moca_blue is currently in development.

Please find more detailed descriptions
of the directories and their role within them, respectively.

This is a pipeline of consecutive operations that can be and will be availabe here.

INPUT DIRECTORY                /0MOTIFS                              /ref_seq
                      - HDF5.file [feature extraction files]       - fastas
                      _____|_________________                      - gffs
                      |                      |                     - meta-data
START DIRECTORY   /mo_nom                   /mo_range                |
output          - get motif patterns      - get motif meta-data      |
                - motif annotation           |                       |
                - motif modification                                 |
                      |______________________________________________|
                      |                                |
                  /mo_clu                             MAPPING to reference (external)
                - analyze motifs             |        use e.g. "blamm
                - compare/cluster            |        (https://github.com/biointec/blamm)
                                             |        cp occurences.txt [results] /mo_proj
                                             |_________|
                                                  |
                                                 /mo_proj
                                                - filter for meaningful matches
                                                - interpret model predictions
                                                - gene annotation
                                                - module generation



mo_nom  --------------------------

Extract motifs from MoDisco hdf5 files and assign nomenclature.
Currently, there are three versions of the same script that can be used for the extraction of a given format of weight matrix.

rdf5_get_xxx_per_pattern.v1.0R.R

PFM - positional frequence matrix
PWM - positional weight matrix (best for clustering/comparison)
CWM - contribution weight matrix (best for mapping)



mo_range ------------------------

Motifs/ EPMs are not distributed at random in a genome.
To optimize the search for motifs/EPMs in a genome or gene-space, these tools
extract the positionally preferred ranges for each motif/EPM in a hdf5 file.

rdf5_get_seql_per_patternV2.R - Extract a list of seqlets and their positions from the hdf5 file

meta_motif_ranges_characteristics_TSS-TTS.1.1.R - Producee a table from the rdf5_get_seql_per_patternV2.R output
  that provides the gene-space statistics for each motif/seqlet in reference to transcription start and stop sites (TSS, TTS)



mo_clu --------------------------

Analyse and Edit motif-files stored in jaspar-format here. Results should be stored in the "out" directory.

mo_cluster_v2.0R - generates dendrograms/trees based on distancy-matrix for different models.


mo_old ---------------------------

Old and outdated scripts used for the moca_blue suite are stored here.



ref_seq -------------------------

Store genome data like fasta, gff and many more here for INPUT. 

