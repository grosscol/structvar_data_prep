#!/bin/bash

REF="/net/share/ftp/vt/grch38/hs38DH.fa"
RANGE="chr1:50185000-50188000"
BAD_FILES=(\
"/net/topmed10/incoming/mapping/results/broad/Meyers/b38/NWD439880/NWD439880.recab.cram" \
"/net/topmed9/working/mapping/results/broad/Gelb/b38/NWD384528/NWD384528.recab.cram")

echo ${BAD_FILES[0]}

samtools view --reference ${REF} ${BAD_FILES[0]} ${RANGE}

