#!/bin/bash
replaceWith=""
for filename in narrowpeak_macs2/*.narrowPeak; do
    sed '1d' "$filename" > "${filename/narrowpeak_macs2\//"$replaceWith"}"
done

