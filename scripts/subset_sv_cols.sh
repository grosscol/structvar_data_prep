#!/usr/bin/env bash

# Expected use
#  ./subset__sv_cols.sh | gzip > struct_var.gz

FMT_JOINED=''
# From: https://stackoverflow.com/a/17841619
function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf -v FMT_JOINED %s "$f" "${@/#/$d}"
  fi
}

# FMT_STR='%CHROM\t%POS\t%REF\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\n'
join_by "\t%" '%CHROM' 'POS' 'REF' 'INFO/SVTYPE' 'INFO/SVLEN' 'INFO/END' 'INFO/AC' 'INFO/AN' 'INFO/PRE' 'INFO/POST'

bcftools query -f "${FMT_JOINED}\n" data/structural.variant.1.1.genotypes.bcf
