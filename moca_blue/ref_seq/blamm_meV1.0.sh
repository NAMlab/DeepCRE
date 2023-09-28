#!/bin/bash

# Define variables
project="Spen_ch01_0e3"
pt_score=0.0001
motifs_file="Soly-MSR_20230118_uu.jaspar"

#
# Start recording runtime and resource information
start_time=$(date +%s)
start_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)


# Run the commands with variables
./blamm dict sequences.mf
./blamm hist -e "$motifs_file" sequences.mf #-e generated empirical PWM scores
./blamm scan -rc -pt "$pt_score" "$motifs_file" sequences.mf

#
echo "CLEANING UP"

mkdir occ$project
mv hist_* occ$project/
mv occurrences.txt occ$project/
mv PWMthresholds.txt occ$project/

#
# Stop recording runtime and resource information
end_time=$(date +%s)
end_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)

# Calculate runtime
runtime=$((end_time - start_time))
# Print runtime and resource information
echo "Script runtime: $runtime seconds"
echo "Resource usage:"
echo "$start_resources" | awk '{print "Start:", $0}'
echo "$end_resources" | awk '{print "End:", $0}'
