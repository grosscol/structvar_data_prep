#!/usr/bin/env bash
#
# Expected use:
#   select_samples.sh
#

SV_BCF=/net/topmed2/working/gt-release/exchange-area/structural.variant.merge/structural.variant.1.1.genotypes.bcf
HETHOM_BIN=/net/wonderland/home/grosscol/projects/het_hom_selector/build/het_hom_selector/bin/het_hom_sel

NO_GENOS_FILE=no_genotypes.tsv
RESULT_FILE=sv_selected.tsv.gz

# Awk script to parse ID column into begin and end -/+ 100 postions
read -r -d '' IDPARSE << "HEREDOC"
BEGIN { FS="\t"; OFS = FS}
/^#/     {
  if($1 == "#CHROM") print $1, $2, $3, $4, $5, "START", "END", $6, $7;
  else print $0;
}
/^chr/   {
  split($3, id_arr,":"); 
  split(id_arr[2], pos_arr, "-");
  bed_start = strtonum(pos_arr[1]) - 100; 
  bed_end = strtonum(pos_arr[2]) + 100;
  
  print "dbg", bed_start, bed_end, "dbg"
  print $1, $2, $3, $4, $5, bed_start, bed_end, $6, $7;
}
HEREDOC


# Awk script to remove entries without genotypes
read -r -d '' NOGENOS << "HEREDOC"
BEGIN { 
 FS="\t"; OFS = FS;
 if(logfile == "") logfile="no_genos_log.tsv"
}
/\t\t$/  {print $0 > logfile} 
!/\t\t$/ {print $0}
HEREDOC

Processing
${HETHOM_BIN} rnd -n 5 --seed 20251009 --emit-id ${SV_BCF}|\
  awk "${IDPARSE}" |\
  awk "${NOGENOS}" |\
  bgzip > "${RESULT_FILE}"
