#!/usr/bin/env bash
#
# Given list of IDs on stdin,
#  Generate tsv with columns: random id, nwd id,  path to cram
#  Optionally take arg for size of random id pool as $1
#
# Expected use:
#   extract_ids.sh | make_mapping.sh

#######################################
# Genereate pool of unique random ids #
#######################################
make_n_rnd_ids(){
  N_IDS=${1:-10}
  LC_ALL=C
  cat /dev/random |\
    tr -dc 'A-Z0-9' |\
    fold -w 10 |\
    awk '!x[$0]++' |\
    head -n ${N_IDS}
}

RND_ID_POOL_SIZE=${1:-150000}
readarray -t RND_IDS < <(make_n_rnd_ids ${RND_ID_POOL_SIZE})

################################
# Map ids to cram path on disk #
################################
MASS_LOOKUP=/net/wonderland/home/grosscol/projects/structvar/pipelines/select/mass_cram_lookup.sh
PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data
RESULTS="${PROC_DATA_DIR}/id_cram_map.tsv"
echo -e "rndid\tnwdid\tpath" > "${RESULTS}"

NLINE=0
BATCH_SIZE=5
INPUT_ARR=()
while read LINE; do
  NLINE=$(( NLINE + 1 ))
  ARR+=($LINE)
  REM=$(( NLINE % BATCH_SIZE ))

  if [ $REM -eq 0 ]; then
    paste <(echo "${RND_IDS[@]:NLINE:BATCH_SIZE}" | tr ' ' $'\n') <("${MASS_LOOKUP}" "${ARR[@]}") >> ${RESULTS}
    ARR=()
  fi
done

if [ $REM -gt 0 ]; then
  paste <(echo "${RND_IDS[@]:NLINE:REM}" | tr ' ' $'\n') <("${MASS_LOOKUP}" "${ARR[@]}") >> ${RESULTS}
fi

NLINE=$(( NLINE + REM))
if [ $NLINE -gt $RND_ID_POOL_SIZE ]; then
  echo "Insufficent random ids generated!" 1>&2
  exit 1
fi

echo "--- done ---"
