#!/bin/bash
replaceWith=""
for filename in *.narrowPeak; do
    python3 processnarrowpeak.py "$filename" > "$filename.processed"
done

