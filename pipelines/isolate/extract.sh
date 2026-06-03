#!/usr/bin/env bash
#
# Extract and compile reads of isolated structural variants.
#
# Read each file of non-overlapping variants and compile a cram of reads
# Input headerless files have tab separated columns: CHROM POS ID REF ALT START END HOM HET


##########
# Inputs #
##########
PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data
INPUT_DIR="${PROC_DATA_DIR}/isolated"

###########
# Outputs #
###########
OUTPUT_DIR=/net/wonderland/home/grosscol/projects/structvar/result/crams
mkdir -p "${OUTPUT_DIR}"

##########
# Config #
##########
MAX_PIPES=100
PIPE_DIR=/tmp/isolate_pipes
MERGE_DIR=/tmp/isolate_merges
mkdir -p ${PIPE_DIR}
mkdir -p ${MERGE_DIR}

#############
# Functions #
#############

generate_output_path(){
  # $1 is output directory.
  # $2 is input index filename.
  local BNAME=$(basename -s .tsv $2)
  echo "${1}/${BNAME}.cram"
}

process_sv_index(){
# $1 Path to a structural variant to input. A crams index file.
  local INFILE=$1
  local N_PIPE=0
  local N_ROUND=0
  local BASE_NAME=$(basename -s '.tsv' "${INFILE}")
  local OUTPUT_CRAM=$( generate_output_path "${OUTPUT_DIR}" "${INFILE}" )

  while read -ra ARR; do
    local CRAM_FILE=${ARR[5]}
    # Replace INV_ DUP_ DEL_ with "chr" in variant id to form region string
    local REGION="chr${ARR[2]:4}"

    # Merge and claer pipes when reach max count of pipes
    if [ $N_PIPE -eq $MAX_PIPES ]; then
      echo "merging pipes to ${MERGE_DIR}/m_${N_ROUND}.cram"
      samtools merge -o "${MERGE_DIR}/m_${N_ROUND}.cram" "${PIPE_DIR}"/*.cram

      rm "${PIPE_DIR}"/*
      sleep 1

      N_PIPE=0
      N_ROUND=$(( N_ROUND + 1 ))
    fi

    # Make new pipe for subsetting associated cram of a structural variant
    local PIPE_NAME="${PIPE_DIR}/pipe_${N_PIPE}.cram"
    mkfifo "${PIPE_NAME}"
    samtools view --cram -o ${PIPE_NAME} ${CRAM_FILE} ${REGION} &

    N_PIPE=$(( N_PIPE + 1 ))
  done < <(head -n 305 $INFILE)

  # Final pipe merge if present
  N_PIPES_REMAIN=$(find "${PIPE_DIR}" -maxdepth 1 -type p | wc -l)

  if [ $N_PIPES_REMAIN -gt 0 ]; then
    samtools merge -o "${MERGE_DIR}/m_${N_ROUND}.cram" ${PIPE_DIR}/*.cram
    rm "${PIPE_DIR}"/*
  fi

  # Merge the merged files into a final cram
  
  local N_MERGE_FILES=$(find "${MERGE_DIR}" -maxdepth 1 -type f | wc -l)

  if [ $N_MERGE_FILES -gt 1 ]; then
    samtools merge -f -o ${OUTPUT_CRAM} ${MERGE_DIR}/*
  elif [ $N_MERGE_FILES -eq 1 ]; then
    local FILES=("${MERGE_DIR}"/*.cram)
    mv ${FILES[0]} ${OUTPUT_CRAM}
  else
    echo "Error: No aggregate crams generated" 1>&2
  fi

  # Clean out merge directory
  rm "${MERGE_DIR}"/*
}

########
# Main #
########

NFILE=0
for FILE in ${INPUT_DIR}/*; do
  # logging
  echo "Processing: ${FILE}"
  DEST=$( generate_output_path "${OUTPUT_DIR}" "${FILE}" )

  # debugging limit run
  NFILE=$(( NFILE + 1 ))
  if [ $NFILE -gt 351 ]; then
    break
  fi

  if stat "${DEST}" > /dev/null 2>&1; then
    echo "Skipping.   (${DEST})"
  else
    process_sv_index $FILE
  fi
done
