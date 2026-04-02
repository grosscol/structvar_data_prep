#!/usr/bin/env bash

# Expected use
#  ./extract_sample_data.sh


# Subset SV calls with genotypes to a couple examples and only GTs that have SV
INFILE="data/structural.variant.1.1.genotypes.bcf"

# Two del examples
# chr1 10309401
# chr11 2624048
bcftools view -r "chr1:10309401" --output-type b --output out.bcf ${INFILE} 
