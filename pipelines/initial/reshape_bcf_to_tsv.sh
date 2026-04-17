#!/usr/bin/env bash


PROC_DATA_DIR=/net/wonderland/home/grosscol/projects/structvar/proc_data
RESULT_DIR=/net/wonderland/home/grosscol/projects/structvar/result

BCF_FULL=/net/topmed2/working/gt-release/exchange-area/structural.variant.merge/structural.variant.1.1.genotypes.bcf
BCF_TAG=${PROC_DATA_DIR}/sv.1.1.genotypes.retagged.bcf

mkdir -p ${PROC_DATA_DIR}

# Re-tag bcf to force AC and AF to match genotypes
#   Remove variants where AC=0
if [ ! -f ${BCF_TAG} ]; then
  bcftools plugin fill-tags ${BCF_FULL} -- -t AN,AC,AC_Het,AC_Hom,AF,F_MISSING,NS |\
    bcftools view --min-ac 1 --output-type b > ${BCF_TAG} 
  bcftools index ${BCF_TAG}
fi


TSV_20=${PROC_DATA_DIR}/sv.1.1.info.chr20.tsv
Q_FMT=\
'%CHROM\t%POS\t%ID\t%REF\t%FILTER\t%INFO/AC\t%INFO/AC_Het\t%INFO/AC_Hom\t%INFO/AF\t%INFO/AN\t%INFO/CALLRATE\t'\
'%INFO/DPCNT\t%INFO/DPCNToverlap\t%INFO/END\t%INFO/F_MISSING\t%INFO/GCPCT\t%INFO/GTCNT\t%INFO/NS\t'\
'%INFO/POST\t%INFO/PRE\t%INFO/SVLEN\t%INFO/SVTYPE\n'

if [ ! -f ${TSV_20} ] || [ ! -s ${TSV_20} ]; then
  bcftools query -r chr20 -f "${Q_FMT}" "${BCF_TAG}" > "${TSV_20}"
fi

echo "done"
