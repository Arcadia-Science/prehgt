#!/bin/bash

files=()

for file in GCF_*; do
    if [[ -f "$file" ]]; then
        key=$(echo "$file" | cut -d '_' -f 2-)
        files+=("$key")
    fi
done

for key in "${files[@]}"; do
    gcf_file="GCF_${key}"
    gca_file="GCA_${key}"
    
    if [[ -f "$gca_file" ]]; then
        rm "$gca_file"
    fi
done
