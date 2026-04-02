#!/usr/bin/env jq

# Trim to range of positions
# Requires command line json arguments: pos, end
# jq -f trimmer.jq --argjson pos 50000 --argjson end 90000 input.json
#
# Doesn't use chr as the dataset is already separated by chromosome

def range_included:
  if type == "object" and .pos != null
  then
    select(.pos > $ARGS.named.pos and .pos < $ARGS.named.end)
  else
    .
  end;

walk(range_included)
