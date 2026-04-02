#!/usr/bin/env bash
#
# Expected use:
#   resolve_samples.sh

INPUT_FILE=sv_selected.tsv.gz
RESULT_FILE=sv_crams.tsv.gz

# where.sh script resovlves cram location from SAMPLE ID
#  requires PROJECT environment variable defined
export PROJECT=topmed
WHERE_SH=/net/mario/cluster/topmed/bin/utils/where.sh

# While reading input file, translate tabs to '@' to avoid consecutive whitespace delimiters.
#  They would be collapsed to a single delimiter due to how Bash does word splitting.
OLD_IFS="${IFS}"
IFS='@'
while read -ra ARR; do
  if [[ "${ARR[0]}" == "#CHROM" ]]; then
    (IFS=$'\t'; echo -e "${ARR[*]}\tCRAMS")
  elif [[ "${ARR[0]}" =~ ^# ]]; then 
    echo -e "${ARR[*]}"
  else 
    IFS=',' read -ra HOMS <<< "${ARR[7]}"
    IFS=',' read -ra HETS <<< "${ARR[8]}"
    IDS=("${HOMS[@]}" "${HETS[@]}")

    ID_STR=$(IFS=','; echo "${IDS[*]}")

    CRAMS=()
    for ID in "${IDS[@]}"; do
      CRAM=$(${WHERE_SH} "${ID}" b38)
      CRAMS+=("${CRAM}")
    done

    CRAMS_STR=$(IFS=';'; echo "${CRAMS[*]}")
    (IFS=$'\t'; echo -e "${ARR[*]}\t${CRAMS_STR}")
  fi
done < <(zcat "${INPUT_FILE}" | tr $'\t' '@') | bgzip > ${RESULT_FILE}
IFS="${OLD_IFS}"
