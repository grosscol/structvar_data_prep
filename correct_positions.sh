#!/usr/bin/env bash

# Correct start and end positions that were expanded in select_samples.sh.
#  The region was expanded by 100bp in each direction to net surrounding
#  relevant reads.
# The correct has been made at the end of the processing pipeline.  This script
#  is to fix the data emitted by the uncorrected pipeline.  The check is easily
#  done by looking if the structural variant id matches the pos and end fields.

zcat rem_input.tgz |\
  jq '. | map(.pos = .pos+100 | .end = .end-100)'
