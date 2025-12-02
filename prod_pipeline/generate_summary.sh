#!/usr/bin/env bash
#
# Process structural variants list with crams annotated.
#   Subsample alignments in variant region for each sample (hets & homs)
#
# Positional params:
#  1. Region pssed to tabix

# Expected use:
#  generate_summary.sh chr1
# Emits summary_sv_reads_chr1.json.gz

REGION="${1:-chr1}"
INPUT_FILE=sv_crams.tsv.gz
CRAM_SUMMARIZER=/net/wonderland/home/grosscol/projects/het_hom_selector/build/cram_summarizer/bin/cram_summ
REF_FILE=/net/topmed/working/atks/ref/hs38DH.fa
OUTPUT_FILE="results/sv_summary_${REGION}.json"

# Results holding directory
mkdir -p pipes
mkdir -p results
echo "[" > ${OUTPUT_FILE}

# While reading input file, translate tabs to '@' to avoid consecutive whitespace delimiters.
#  They would be collapsed to a single delimiter due to how Bash does word splitting.
OLD_IFS="${IFS}"
IFS='@'
ROWNUM=0
while read -ra ROW; do
  if [[ ${ROWNUM} -gt 0 ]]; then
    echo ',' >> ${OUTPUT_FILE}
  fi

  IFS=';' read -ra CRAMS <<< "${ROW[9]}"
  CHR=${ROW[0]}
  START=${ROW[5]}
  END=${ROW[6]}
  SV_ID=${ROW[2]}
  DOWNSAMPLE=$(( 5 * ${#CRAMS[@]} ))

  echo "Processing ${ROWNUM} C:${#CRAMS[@]} S:${DOWNSAMPLE}"

  # Generate one pipe per cram and begin reading into it.
  PIPES=()
  PIPE_COUNT=0
  for CRAM in "${CRAMS[@]}"; do
    if [ -f "${CRAM}" ]; then
      PIPE_NAME="pipes/pipe_${REGION}_${PIPE_COUNT}.bam"
      mkfifo ${PIPE_NAME}
      PIPES+=( ${PIPE_NAME} )
      let "PIPE_COUNT = PIPE_COUNT + 1"

      samtools view -x XA -b --reference ${REF_FILE} ${CRAM} "${CHR}:${START}-${END}" > ${PIPE_NAME} &
    else
      echo "CRAM not found! ${CRAM}" >&2
    fi
  done

  # # Merge all crams into a single pipe
  PIPE_MERGE="pipes/pipe_${REGION}_merge.bam"
  mkfifo "${PIPE_MERGE}"
  samtools merge -f -o "${PIPE_MERGE}" "${PIPES[@]}" &

  # Summarize merged bams, remove unusable reads, and annote.
  # Reduce data size by reducing key name and using 1|0 instead of boolean
  # Append to output file
  ${CRAM_SUMMARIZER} -s "${DOWNSAMPLE}" "${PIPE_MERGE}" |\
    jq --compact-output --arg sv_id "${SV_ID}" --arg chr "${CHR}" \
      --argjson pos ${START} --argjson stop ${END} -f trim_by_chrom.jq |\
    sed 's/is_reverse/rev/;s/true/1/;s/false/0/' >> ${OUTPUT_FILE}

  # Cleanup pipes
  for P in ${PIPES[@]}
  do
    rm ${P}
  done
  rm "${PIPE_MERGE}"

  let "ROWNUM = ROWNUM + 1"
done < <(tabix "${INPUT_FILE}" "${REGION}" | tr $'\t' '@')

echo -n "]" >> ${OUTPUT_FILE}

echo "gzipping results to ${OUTPUT_FILE}.gz"
gzip -f ${OUTPUT_FILE}
