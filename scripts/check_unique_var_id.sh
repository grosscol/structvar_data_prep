#!/usr/bin/env bash

# Check if chr-pos-ref-alt (cpra) is sufficient for identifying structvars.
# Facilitate looking for duplicates

# Expected use:
#    ./check_unique_var_id.sh input.bcf > variant_id_unique_counts.txt


bcftools query -f '%ID%POS%REF%ALT\n' ${1} | sort | uniq -c
