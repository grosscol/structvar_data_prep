#!/usr/bin/env bash

# Experimental data processing
#
# Expected use
#  ./count_het_homs.sh path/to/variants.tgz > het_hom_counts.txt


if [ "$#" -ne 1 ]
then
  echo "Need path to variants.tgz"
  exit 1
fi

# decompress
# cut to het & hom columns
# delete anything that's not a comma or tab (field separator)
# skip first four rows then count the length of each field (homs & hets)
TAB=$'\t'
zcat ${1} |\
 cut -f 6-7 |\
 sed "s/[^,${TAB}]//g" |\
 awk -F $'\t' 'NR >4 {print length($1)+1, length($2)+1}'
