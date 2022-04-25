#!/bin/bash
replaceWith=""
for filename in fasta_files/*.gz; do
    python3 negativetrain.py "$filename" 
done

