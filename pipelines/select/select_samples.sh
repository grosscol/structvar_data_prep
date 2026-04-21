#!/usr/bin/env bash
#
# Expected use:
#   select_samples.sh

# Source Data
PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data
BCF_TAG=${PROC_DATA_DIR}/sv.1.1.genotypes.retagged.bcf

# Tooling
HETHOM_BIN=/net/wonderland/home/grosscol/projects/het_hom_selector/build/het_hom_selector/bin/het_hom_sel

# In process intermediate output
NO_GENOS_FILE=${PROC_DATA_DIR}/no_genotypes.tsv
SELECTED_IDS_FILE=${PROC_DATA_DIR}/sv_selected.tsv.gz

# Awk script to parse ID column into bed_start and bed_end -/+ 100 postions
read -r -d '' IDPARSE << "HEREDOC"
BEGIN {FS="\t"; OFS = FS}
/^#/     {
  if($1 == "#CHROM") print $1, $2, $3, $4, $5, "START", "END", $6, $7;
  else print $0;
}
/^chr/   {
  split($3, id_arr,":");
  split(id_arr[2], pos_arr, "-");
  bed_start = strtonum(pos_arr[1]) - 100;
  bed_end = strtonum(pos_arr[2]) + 100;

  print $1, $2, $3, $4, $5, bed_start, bed_end, $6, $7;
}
HEREDOC


# Awk script to remove entries without genotypes
#  Identified by two consecutive delimters at end of file indicating:
#  no het or hom genotypes called.
# This file should be empty as the initial processing should have removed these variants.
read -r -d '' NOGENOS << "HEREDOC"
BEGIN {
 FS="\t"; OFS = FS;
 if(logfile == "") logfile="no_genos_log.tsv"
}
/\t\t$/  {print $0 > logfile}
!/\t\t$/ {print $0}
HEREDOC

# Processing
${HETHOM_BIN} rnd -n 1 --seed 20251009 --emit-id ${BCF_TAG}|\
  awk "${IDPARSE}" |\
  awk "${NOGENOS}" |\
  bgzip > "${SELECTED_IDS_FILE}"
