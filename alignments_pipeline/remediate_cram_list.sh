#!/usr/bin/env bash

# This script corrects the data with the paths to crams as the final column.
# The error was in handling rows where HET or HOM was blank and the two cols were collapsed.

SELECTED_TSV=sv_selected.tsv.gz
CRAMS_LIST=crams.txt
OUTFILE=out.tsv.gz

read -r -d '' HEADER << "HEREDOC"
#RANDOM_SEED=20251009
#MAX_RANDOM_HOM_HETS=<<5
#SAMPLES_USED=NA
#CHROM\tPOS\tID\tREF\tALT\tSTART\tEND\tHOM\tHET\tCRAMS
HEREDOC

# Input pipes
PIPE1=pipe_one
PIPE2=pipe_two
mkfifo ${PIPE1}
mkfifo ${PIPE2}

# Setup output pipe for writing
PIPEO=pipe_out
mkfifo ${PIPEO}
bgzip < ${PIPEO} > $OUTFILE &

# Start writing to input pipes
(zcat $SELECTED_TSV | tail -n +5 > $PIPE1) &
(tail -n +5 $CRAMS_LIST > $PIPE2) &

# Add header and combine inputs into output pipe
exec 3>${PIPEO}
echo -e "${HEADER}" > ${PIPEO}
paste ${PIPE1} ${PIPE2} | column -s $'\t' -o $'\t' > ${PIPEO}
exec 3>&-

rm ${PIPE1}
rm ${PIPE2}
rm ${PIPEO}
