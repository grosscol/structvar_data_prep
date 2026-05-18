#!/usr/bin/env bash
#
# Make lists of non-overlapping structural variants with at least 50 bp between them.
# Read throught the file, splitting input lines into separate output files when the struct variant
# would overlap with a previous structural variant.

IFS=$'\t'
PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data

##########
# Inputs #
##########
INPUT_FILE=${PROC_DATA_DIR}/sv_selected.tsv.gz

###########
# Outputs #
###########
OUTPUT_DIR="${PROC_DATA_DIR}/isolated"
INDEX_FILE=${PROC_DATA_DIR}/sv_cram_map.tsv

mkdir -p ${OUTPUT_DIR}

###############################################
# Derermine where header ends and data begins #
###############################################
DATA_START_ROW=1
while read LINE; do
  if [[ ! "${LINE}" =~ ^# ]]; then
    break
  fi
  DATA_START_ROW=$(( DATA_START_ROW + 1 ))
done < <(zcat "${INPUT_FILE}")

#############################################
# State: ~Tracking output buffers and files #
#############################################
# File destinations, end positions, and buffers for output data
DEST_ARR=()
ENDS_ARR=()
BUFF_ARR=()
# Min number of positions (base pairs) to leave bewteen structvars
BP_SPACE=200

# Current chromosome and index into destination array and buffer
CURR_CHR=''
CURR_DEST_IDX=''
RECORD_NUM=0

#####################
# Support functions #
#####################
generate_new_dest(){
  # Args: DIR, CHR, INDEX
  printf -v DEST_STR "%s_selected_p%.3d.tsv" $2 $3
  echo "${1}/${DEST_STR}"
}

flush_dest_buffers(){
  for IDX in "${!DEST_ARR[@]}"; do
    echo -e "${BUFF_ARR[$IDX]}" > "${DEST_ARR[$IDX]}"
  done
  # Clear arrays
  DEST_ARR=()
  ENDS_ARR=()
  BUFF_ARR=()
}

debug_state(){
  echo "DEST_ARR length ${#DEST_ARR[@]}"
  echo "ENDS_ARR length ${#ENDS_ARR[@]}"
  echo "BUFF_ARR length ${#BUFF_ARR[@]}"
  echo "RECORD_NUM: $RECORD_NUM"
  echo "$CHR $START $END $ID"
}

#############
# Main loop #
#############
while read -ra ARR; do
  CHR=${ARR[0]}
  START=${ARR[5]}
  END=${ARR[6]}
  ID=${ARR[2]}

  # Clear the buffers and start over when new chromosome is encountered
  if [ ! "${CHR}" = "${CURR_CHR}" ]; then
    flush_dest_buffers

    CURR_CHR=${CHR}
    CURR_DEST_IDX=0

    DEST_ARR[${CURR_DEST_IDX}]=$(generate_new_dest "${OUTPUT_DIR}" "${CHR}" "${CURR_DEST_IDX}")
    echo "" > ${DEST_ARR[${CURR_DEST_IDX}]}

    ENDS_ARR[$CURR_DEST_IDX]=0
    BUFF_ARR[$CURR_DEST_IDX]=""
  fi

  # When overlap with current buffer and start position of sv
  # Find existing buffer withouth overlap or begin a new buffer
  if [ ! ${START} -gt ${ENDS_ARR[$CURR_DEST_IDX]} ]; then
    LEN=${#ENDS_ARR[@]}
    PREV_IDX=$((CURR_IDX - 1))
    NEXT_IDX=$((CURR_IDX + 1))
    END_IDX=$(( LEN - 1))

    # Search forward through buffers
    LOWER_IDXS=$(seq -s ' '  0 ${PREV_IDX})
    UPPER_IDXS=$(seq -s ' '  ${NEXT_IDX} ${END_IDX})
    SEARCH_ORDER=(${UPPER_IDXS} ${LOWER_IDXS})

    CURR_DEST_IDX=-1
    for IDX in "${!ENDS_ARR[@]}"; do
      if [ ${START} -gt ${ENDS_ARR[$IDX]} ]; then
        CURR_DEST_IDX=$IDX
        break
      fi
    done

    # Create new buffer when non-overlapping is not available
    if [ $CURR_DEST_IDX = -1 ]; then
      LEN=${#ENDS_ARR[@]}
      CURR_DEST_IDX=${LEN}

      ENDS_ARR[${CURR_DEST_IDX}]=0
      DEST_ARR[${CURR_DEST_IDX}]=$(generate_new_dest "${OUTPUT_DIR}" "${CHR}" "${CURR_DEST_IDX}")
      BUFF_ARR[${CURR_DEST_IDX}]=""
      CURR_DEST_IDX=$NEXT_DEST_IDX
    fi
  fi

  # Append line to output buffer.
  BUFF_ARR[${CURR_DEST_IDX}]+="${ARR[*]}\n"
  ENDS_ARR[${CURR_DEST_IDX}]=$((END + BP_SPACE))

  # Emit corresponding index record
  echo -e "${ID}\t${CURR_DEST_IDX}" >> "${INDEX_FILE}"

  RECORD_NUM=$(( RECORD_NUM + 1 ))
done < <(zcat "${INPUT_FILE}" | tail -n +${DATA_START_ROW} | head -n 500)

# Final flush of buffers to disk
flush_dest_buffers
