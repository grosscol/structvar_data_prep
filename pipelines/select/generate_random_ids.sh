#!/usr/bin/env bash
#
# Generate list of unique ids to obfuscate actual sample ids.
# Use:
#   generate_random_ids.sh N_IDS ID_LEN > rnd_ids_file.txt

# Default values for inputs
N_IDS=${1:-200}
ID_LEN=${2:-12}

# Validate inputs are numeric
NUM_PATT='[0-9]+$'
if ! [[ $N_IDS =~ $NUM_PATT ]]; then
  echo "N_IDS needs to be positive integer" 1>&2;
  exit 1
fi

if ! [[ $ID_LEN =~ $NUM_PATT ]]; then
  echo "ID_LEN needs to be positive integer" 1>&2;
  exit 1
fi

# Compute variables
let "N_RND = N_IDS + 100"
let "N_RAW_CHARS = N_RND * ID_LEN"

#readarray -t RND_IDS < <( cat /dev/urandom |\
#  LC_ALL=C tr -dc A-Z0-9 |\
#  head -c ${N_RAW_CHARS} |\
#  fold -w ${ID_LEN} - |\
#  awk '!x[$0]++')

cat /dev/urandom |\
  LC_ALL=C tr -dc A-Z0-9 |\
  head -c ${N_RAW_CHARS} |\
  fold -w ${ID_LEN} - |\
  awk '!x[$0]++'
