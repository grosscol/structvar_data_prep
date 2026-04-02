#!/usr/bin/env bash

# Check if chr-pos-ref-alt (cpra) is sufficient for identifying structvars.
# Facilitate looking for duplicates

# Expected use:
#    ./check_unique_cpra_id.sh input.bcf > cpra_id_unique_counts.txt


bcftools view ${1} | grep -v '^#' | cut -f 1,2,4,5 | tr -d '\t' | sort | uniq -c
