

project="Spen_ch01_0e2"
inputFile="occurrences.txt"
outputSize="1000000"

# Start recording runtime and resource information
start_time=$(date +%s)
start_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)

split -l $outputSize --numeric-suffixes occ$project/$inputFile smallfile
#mkdir occ$project
mv smallfile* occ$project/
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
