#!/bin/bash
replaceWith=""
for filename in noheader_narrowpeaks/*.narrowPeak; do
    python3 processnarrowpeak.py "$filename" > "$filename.stats"
done

