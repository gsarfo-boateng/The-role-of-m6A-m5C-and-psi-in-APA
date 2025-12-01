#!/bin/bash

for file in *_PAS_data.txt; do
    bed_file="${file%.txt}.bed"
    awk 'NR > 1 { print $1, $2, $3 }' OFS="\t" "$file" > "$bed_file"
    echo "Created $bed_file"
done
