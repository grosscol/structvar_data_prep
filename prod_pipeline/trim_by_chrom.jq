#!/usr/bin/env jq

# Command line passed in variables: chrom, pos, end, id 


def chroms_only:
  if type == "object" and .chr != null
  then 
    select(.chr | test("^chr([0-9]{1,2}|X|Y)$"))
  else 
    . 
  end;

walk(chroms_only) | .sv_id=$sv_id | .chr=$chr | .pos=$pos | .end=$stop
