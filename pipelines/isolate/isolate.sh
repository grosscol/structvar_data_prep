#!/usr/bin/env bash
#
# Make lists of non-overlapping structural variants with at least 50 bp between them.
# Read throught the file, splitting input lines into separate output files when the struct variant
# would overlap with a previous structural variant.
#
# Output headerless files with tab separated columns: CHROM POS ID REF ALT START END HOM HET

IFS=$'\t'
PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data

##########
# Inputs #
##########
INPUT_FILE=${PROC_DATA_DIR}/sv_selected.tsv.gz
ID_CRAM_FILE=${PROC_DATA_DIR}/id_cram_map.tsv

###########
# Outputs #
###########
OUTPUT_DIR="${PROC_DATA_DIR}/isolated"
INDEX_FILE=${PROC_DATA_DIR}/sv_cram_map.tsv
mkdir -p ${OUTPUT_DIR}
truncate -s 0 ${INDEX_FILE}

###################################
# Read id -> cram map into memory #
###################################

declare -A ID_CRAM_ARR
while read -ra ARR; do
   ID=${ARR[1]}
   VAL=${ARR[2]}
   ID_CRAM_ARR[$ID]=${VAL}
done < <(tail -n +2 ${ID_CRAM_FILE})

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
# File destinations, end positions, output records, and index of what went where
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
  # Args: CHR, INDEX
  printf "%s_selected_p%.3d.tsv" $1 $2
}

flush_dest_buffers(){
  for IDX in "${!DEST_ARR[@]}"; do
    # Chomp last two chracters "\n" to avoid writing trailing endline at end of file
    echo -e "${BUFF_ARR[$IDX]::-2}" > "${OUTPUT_DIR}/${DEST_ARR[$IDX]}"
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

    DEST_FILE=$(generate_new_dest "${CHR}" "${CURR_DEST_IDX}")
    DEST_ARR[${CURR_DEST_IDX}]="${DEST_FILE}"
    echo "" > "${OUTPUT_DIR}/${DEST_ARR[${CURR_DEST_IDX}]}"

    ENDS_ARR[$CURR_DEST_IDX]=0
    BUFF_ARR[$CURR_DEST_IDX]=""
  fi

  # When overlap with current buffer and start position of sv
  # Find existing buffer withouth overlap or begin a new buffer
  if [ ! ${START} -gt ${ENDS_ARR[$CURR_DEST_IDX]} ]; then
    LEN=${#ENDS_ARR[@]}
    PREV_IDX=$((CURR_DEST_IDX - 1))
    NEXT_IDX=$((CURR_DEST_IDX + 1))
    END_IDX=$(( LEN - 1))

    # Search forward through buffers
    LOWER_IDXS=$(seq -s "$IFS"  0 ${PREV_IDX})
    UPPER_IDXS=$(seq -s "$IFS"  ${NEXT_IDX} ${END_IDX})
    SEARCH_ORDER=(${UPPER_IDXS} ${LOWER_IDXS})

    CURR_DEST_IDX=-1
    for IDX in "${SEARCH_ORDER[@]}"; do
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
      DEST_FILE=$(generate_new_dest "${CHR}" "${CURR_DEST_IDX}")
      DEST_ARR[${CURR_DEST_IDX}]="${DEST_FILE}"
      BUFF_ARR[${CURR_DEST_IDX}]=""
    fi
  fi

  # Lookup CRAM FILE by ID.  Prefer id of homozygous sample.
  if [ -n "${ARR[7]}" ]; then
    SAMPLE_ID=${ARR[7]}
  else
    SAMPLE_ID=${ARR[8]}
  fi
  CRAM_FILE=${ID_CRAM_ARR[$SAMPLE_ID]}

  # Generate line and append to output buffer: CHR POS ID START STOP CRAMFILE
  BASE_INFO="${ARR[@]:0:3}\t${ARR[@]:5:2}"
  OUT_LINE="${BASE_INFO}\t${CRAM_FILE}\n"
  BUFF_ARR[${CURR_DEST_IDX}]+="$OUT_LINE"

  # Update buffer's tracked end position.
  ENDS_ARR[${CURR_DEST_IDX}]=$((END + BP_SPACE))

  # Emit corresponding index record mapping structvar ids to cram file sequences will be in.
  echo -e "${BASE_INFO}\t${DEST_ARR[${CURR_DEST_IDX}]}" >> "${INDEX_FILE}"

  RECORD_NUM=$(( RECORD_NUM + 1 ))
done < <(zcat "${INPUT_FILE}" | tail -n +${DATA_START_ROW})

# Final flush of buffers to disk
flush_dest_buffers
