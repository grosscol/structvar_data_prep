#!/usr/bin/env bash

# Expected use
#  ./extract_het_homs.sh
#
# Creates structvar_het_homs.tgz in working directory.
# Uses the het_hom_selector program written for this purpose

INFILE="data/structural.variant.1.1.genotypes.bcf"
OUTFILE="structvar_het_homs.tgz"

het_hom_selector all --emit-ids --file ${INFILE} |\
  awk '!/\t\t$' |\
  bgzip > ${OUTFILE}
