#!/usr/bin/env bash

# Get a list of samples from structural variants.

# bcftools view -H -r chr11:9207515 data/structural.variant.1.1.genotypes.bcf | head
REGIONS=chr11:2537753
# bcftools query -f '[%SAMPLE %GT\n]' -r ${REGIONS} data/structural.variant.1.1.genotypes.bcf |\
#   cut -d ' ' -f 2 | sort | uniq -c

bcftools view -H -r $REGIONS data/structural.variant.1.1.genotypes.bcf | head
# {
#   "ac": 56,
#   "an": 270994,
#   "chrom": "11",
#   "end": 10297271,
#   "pos": 2537753,
#   "post": "1.07,0.082615",
#   "pre": "1.02,0.084702",
#   "ref": "G",
#   "sv_len": 7759519,
#   "sv_type": "DUP"
# },
# {
#   "ac": 1,
#   "an": 275502,
#   "chrom": "11",
#   "end": 9207515,
#   "pos": 2693282,
#   "post": "1.03,0.080546",
#   "pre": "1.05,0.078697",
#   "ref": "C",
#   "sv_len": 6514234,
#   "sv_type": "DEL"
# },
# {
#   "ac": 1,
#   "an": 276268,
#   "chrom": "11",
#   "end": 5387994,
#   "pos": 5016392,
#   "post": "1.00,0.000000",
#   "pre": "1.00,0.000000",
#   "ref": "T",
#   "sv_len": 371603,
#   "sv_type": "INV"
# },
