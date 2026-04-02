#!/usr/bin/env bash
#
# Trim structvar alignment summary to a specified region

# --compact-output
INFILE=sv_summary_chr11.json.gz
zcat "${INFILE}" | jq -f trim_by_region.jq > /tmp/cgout.json
