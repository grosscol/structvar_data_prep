#!/usr/bin/env bash
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
process_sv_index(){
# $1 Path to a structural variant to crams index file
  local N_PIPE=0
  local N_ROUND=0
  local BASE_NAME=$(basename -s '.tsv' $1)

  while read -ra ARR; do
    local CRAM_FILE=${ARR[5]}
    local REGION="${ARR[0]}:${ARR[2]:6}"

    # Merge and claer pipes when reach max count of pipes
    if [ $N_PIPE -eq $MAX_PIPES ]; then
      echo "Merge and clear the pipes!"
      samtools merge -o "${MERGE_DIR}/m_${N_ROUND}.cram" "${PIPE_DIR}"/*.cram
      rm "${PIPE_DIR}"/*
      sleep 1
      N_PIPE=0
    fi

    # Make new pipe for subsetting associated cram of a structural variant
    local PIPE_NAME="${PIPE_DIR}/pipe_${N_PIPE}.cram"
    mkfifo "${PIPE_NAME}"
    echo "filling pipe $PIPE_NAME"
    samtools view --cram -o ${PIPE_NAME} ${CRAM_FILE} ${REGION} &

    N_PIPE=$(( N_PIPE + 1 ))
    N_ROUND=$(( N_ROUND + 1 ))
  done < <(head -n 305 $1)

  # Final pipe merge if present
  N_PIPES_REMAIN=$(find "${PIPE_DIR}" -maxdepth 1 -type p | wc -l)
  echo "Pipes remaining: $N_PIPES_REMAIN"

  if [ $N_PIPES_REMAIN -gt 0 ]; then
    echo "Merging last ${N_PIPES_REMAIN} pipes"
    samtools merge -o "${MERGE_DIR}/m_${N_ROUND}.cram" ${PIPE_DIR}/*.cram
    rm "${PIPE_DIR}"/*
  fi

  # Merge the merged files into a final cram
  local OUTPUT_CRAM="${OUTPUT_DIR}/${BASE_NAME}.cram"
  local N_MERGE_FILES=$(find "${MERGE_DIR}" -maxdepth 1 -type f | wc -l)
  if [ $N_MERGE_FILES -gt 1 ]; then
    echo "Final merge of subset crams"
    samtools merge -f -o ${OUTPUT_CRAM} ${MERGE_DIR}/*
  elif [ $N_MERGE_FILES -eq 1 ]; then
    echo "One cram generated.  No merge required. Moving to results"
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

# NFILE=0
# for FILE in ${INPUT_DIR}/*; do
#   BASE_NAME=$(basename -s '.tsv' $FILE)
#   NFILE=$(( NFILE + 1 ))
#   if [ $NFILE -gt 5 ]; then
#     break
#   fi
#   process_sv_index $FILE
# done

# Biggest file
# INFILE=/net/wonderland/home/grosscol/projects/structvar/proc_data/isolated/chr2_selected_p037.tsv
INFILE=/net/wonderland/home/grosscol/projects/structvar/proc_data/isolated/chr9_selected_p000.tsv
echo "debug"
date
wc -l $INFILE
process_sv_index $INFILE
date

# debug
# NFILE=/net/topmed3/working/mapping/results/washu/BioMe-COPD/b38/NWD411393/NWD411393.recab.cram
# REGION="chr1:1098501-1110400"
# 
# mkfifo /tmp/cg_fifo.sam
# samtools view -o /tmp/out.cram --cram /tmp/cg_fifo.sam &
# samtools view --cram "${NFILE}" "${REGION}" > /tmp/cg_fifo.cram
# samtools view --cram "${NFILE}" "${REGION}" >> /tmp/cg_fifo.cram
# 
# for IDX in $(seq 1 3); do
# done

