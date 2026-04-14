#!/usr/bin/env bash

BCF_FULL=/net/topmed2/working/gt-release/exchange-area/structural.variant.merge/structural.variant.1.1.genotypes.bcf
SCRATCH=/net/wonderland/home/grosscol/projects/structvar/scratch

BCF_20=${SCRATCH}/structural.variant.chr20.bcf
BCF_20_TAG=${SCRATCH}/structural.variant.chr20.retag.bcf
BCF_20_TAG_UNC=${SCRATCH}/structural.variant.chr20.retag.uncalled.bcf
BCF_20_TAG_AC1=${SCRATCH}/structural.variant.chr20.retag.minac1.bcf
BCF_20_TAG_AC0=${SCRATCH}/structural.variant.chr20.retag.maxac0.bcf

# Subset to chr to if not present

if [ ! -f ${BCF_20} ]; then
  bcftools view -r chr20 --output-type b ${BCF_FULL} > ${BCF_20}
  bcftools index ${BCF_20}
fi

if [ ! -f ${BCF_20_TAG} ]; then
  bcftools plugin fill-tags --output-type b ${BCF_20} -- -t AN,AC,AC_Het,AC_Hom,AF,F_MISSING,NS > ${BCF_20_TAG} 
  bcftools index ${BCF_20_TAG}
fi

if [ ! -f ${BCF_20_TAG_UNC} ]; then
  bcftools view -r chr20 --output-type b --uncalled ${BCF_20_TAG} > ${BCF_20_TAG_UNC}
  bcftools index ${BCF_20_TAG_UNC}
fi

if [ ! -f ${BCF_20_TAG_AC1} ]; then
  bcftools view -r chr20 --output-type b --min-ac 1 ${BCF_20_TAG} > ${BCF_20_TAG_AC1}
  bcftools index ${BCF_20_TAG_AC1}
fi

if [ ! -f ${BCF_20_TAG_AC0} ]; then
  bcftools view -r chr20 --output-type b --max-ac 0 ${BCF_20_TAG} > ${BCF_20_TAG_AC0}
  bcftools index ${BCF_20_TAG_AC0}
fi
stat ${BCF_20_TAG_AC1}
bcftools view -H ${BCF_20_TAG_AC1} | wc -l
bcftools view -H ${BCF_20_TAG_AC0} | wc -l

echo "done"
