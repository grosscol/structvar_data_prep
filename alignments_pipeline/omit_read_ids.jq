#!/usr/bin/env jq

def chroms_only:
  if type == "object" and .chr != null
  then 
    select(.chr | test("^chr([0-9]{1,2}|X|Y)$"))
  else 
    . 
  end;

walk(chroms_only) | .sv_id="foo" | .chr="chrB" | .pos=200 | .end=400 |
   .all_splits=(.all_splits | map(.)) | .all_pairs=(.all_pairs | map(.)) 
