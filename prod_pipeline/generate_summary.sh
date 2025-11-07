#!/usr/bin/env bash

INPUT_FILE=sv_crams.tsv.gz
CRAM_SUMMARIZER=/net/wonderland/home/grosscol/projects/het_hom_selector/build/cram_summarizer/bin/cram_summ
REF_FILE=/net/topmed/working/atks/ref/hs38DH.fa

# While reading input file, translate tabs to '@' to avoid consecutive whitespace delimiters.
#  They would be collapsed to a single delimiter due to how Bash does word splitting.
OLD_IFS="${IFS}"
IFS='@'
SAMPLE=0
while read -ra ROW; do
  echo "--------------------- ${SAMPLE}"
  IFS=';' read -ra CRAMS <<< "${ROW[9]}"
  CHR=${ROW[0]}
  START=${ROW[5]}
  END=${ROW[6]}

  # Generate one pipe per cram and begin reading into it.
  PIPES=()
  PIPE_COUNT=0
  for CRAM in "${CRAMS[@]}"; do
    if [ -f "${CRAM}" ]; then
      PIPE_NAME="pipe_${PIPE_COUNT}.bam"
      mkfifo ${PIPE_NAME}
      PIPES+=( ${PIPE_NAME} )
      let "PIPE_COUNT = PIPE_COUNT + 1"

      samtools view -x XA -b --reference ${REF_FILE} ${CRAM} "${CHR}:${START}-${END}" > ${PIPE_NAME} &
    else
      echo "CRAM not found! ${CRAM}" >&2
    fi

  done

  mkfifo pipe_out.bam

  # # Merge all crams
  # #samtools merge -f -o pipe_out.bam "${PIPES[@]}" &
  echo "Filling merge pipe"
  # samtools merge -f -o "out_${SAMPLE}.bam" "${PIPES[@]}" 
  samtools merge -f -o "pipe_out.bam" "${PIPES[@]}" &

  # # Summarize merged bams.
  echo "Summarizing"
  # ${CRAM_SUMMARIZER} "out_${SAMPLE}.bam" | gzip > "summary_${SAMPLE}.json.gz"
  ${CRAM_SUMMARIZER} "pipe_out.bam" | gzip > "summary_${SAMPLE}.json.gz"

  # Cleanup pipes
  for P in ${PIPES[@]}
  do
    rm ${P}
  done
  rm pipe_out.bam

  let "SAMPLE = SAMPLE + 1"
done < <(zcat "${INPUT_FILE}" | tail -n +5 | head -n 5 | tr $'\t' '@')
