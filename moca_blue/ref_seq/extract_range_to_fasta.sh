#!/bin/bash

#It extracts the gene ranges from the GFF file and saves them to the extracted_ranges.txt file.
#After that, it uses samtools to extract the corresponding sequences from the FASTA file (fasname) based on the ranges provided in output_file.
#The extracted sequences are then saved to fasout.

#20230628 Simon M. Zumkeller
# Specify the integer value to subtract/add from feature start and end
# ALL features are extracted from the leading strand
flank_size=1000
####################################################
filename="ITAG4.0_gene_models.gff"  # Replace with the path to your GFF file
fasname="S_lycopersicum_chromosomes.4.00.fa" # Replace with the path to your fas file
####################################################
output_file="extracted_ranges.txt"
#####################################################
fasout="${fasname%.*}_1kbp-flank.fa"



while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9; do
    if [[ $col3 == *"gene"* ]]; then
        col4=$((col4 - flank_size))
        col5=$((col5 + flank_size))
        echo -e "$col1:$col4-$col5"
    fi
done < "$filename" > "$output_file"

samtools faidx "$fasname" -o "$fasout" -r "$output_file" --mark-strand rc
