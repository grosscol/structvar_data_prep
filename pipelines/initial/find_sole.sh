
SCRATCH=/net/wonderland/home/grosscol/projects/structvar/scratch
BCF_20_TAG=${SCRATCH}/structural.variant.chr20.retag.bcf

POS=( 11493816 11822378 12075990 12076001 12190620 12190763 12208901 12231906 12253501 12464771 12537866 \
      13500183 13652772 13699353 13734004 14158489 14264801 14288310 14354296 14358182 14695596 14722701 \
      14728772 14791916 14821682 14852093 14943045 14957671 15030857 15205794 15378201 16018101 16115401 \
      16241237 16325201 16537502 17786001 18134228 18192508 19573050 20602840 20976141 21495870 21608157 \
      21833727 22709771 23357076 23444004 )

# for P in ${POS[*]}; do
#   echo -n "${P}: "
#   bcftools view -H -r "chr20:${P}" ${BCF_20_TAG} | wc -l
# done

# Manual check counts of called genotype
P=${POS[1]}
ID="DUP_20:11822378-12307516"
echo $P
bcftools view -H -r "chr20:${P}" ${BCF_20_TAG} | grep --fixed-strings "${ID}" > /tmp/sv_line.txt
grep -o --fixed-strings -e "./." -e "0/0" -e "1/0" -e "0/1" -e "1/1" /tmp/sv_line.txt | sort | cut -d : -f 1 | uniq -c


bcftools view --output-type b -r "chr20:${P}" ${BCF_20_TAG} > /tmp/sv_line.bcf
bcftools index /tmp/sv_line.bcf
bcftools view -H --uncalled /tmp/sv_line.bcf
