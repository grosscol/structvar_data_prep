#!/usr/bin/env bash
# Extract, sort, and uniq the selected ids
#
# Expected use:
#   extract_ids.sh > extracted_ids.txt

PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data
INPUT_FILE=${PROC_DATA_DIR}/sv_selected.tsv.gz

# Derermine where header ends and data begins
DATA_START_ROW=1
while read LINE; do
  if [[ ! "${LINE}" =~ ^# ]]; then
    break
  fi
  DATA_START_ROW=$(( DATA_START_ROW + 1 ))
done < <(zcat "${INPUT_FILE}")

# Read file starting with the data rows.
#  extract the ID columns. combine tab separated cols into one column.
#  remove blank lines and sort ids
zcat ${INPUT_FILE} | tail -n +${DATA_START_ROW} |\
  cut -f 8,9 | tr $'\t' $'\n' |\
  awk NF | sort | uniq
