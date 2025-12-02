#!/usr/bin/env jq

# Command line passed in variables: chrom, pos, end, id
# Reduce file size by:
#   Remove reads mapped to non-cannonical chromosomes
#   Omit alignment ids in favor of flat array of alignments

def chroms_only:
  if type == "object" and .chr != null
  then
    select(.chr | test("^chr([0-9]{1,2}|X|Y)$"))
  else
    .
  end;

walk(chroms_only) | .sv_id=$sv_id | .chr=$chr | .pos=$pos | .end=$stop |
  .all_splits=(.all_splits | map(.)) | .all_pairs=(.all_pairs | map(.))
