#!/usr/bin/env bash
nohup ./generate_summary.sh chr1  > results/chr1.log &
nohup ./generate_summary.sh chr2  > results/chr2.log &
nohup ./generate_summary.sh chr3  > results/chr3.log &
nohup ./generate_summary.sh chr4  > results/chr4.log &
nohup ./generate_summary.sh chr5  > results/chr5.log &

# nohup ./generate_summary.sh chr6  > results/chr6.log &
# nohup ./generate_summary.sh chr7  > results/chr7.log &
# nohup ./generate_summary.sh chr8  > results/chr8.log &
# nohup ./generate_summary.sh chr9  > results/chr9.log &
# nohup ./generate_summary.sh chr10 > results/chr10.log &
# nohup ./generate_summary.sh chr11 > results/chr11.log &

#nohup ./generate_summary.sh chr12 > results/chr12.log &
#nohup ./generate_summary.sh chr13 > results/chr13.log &
#nohup ./generate_summary.sh chr14 > results/chr14.log &
#nohup ./generate_summary.sh chr15 > results/chr15.log &
#nohup ./generate_summary.sh chr16 > results/chr16.log &
 
# nohup ./generate_summary.sh chr17 > results/chr17.log &
# nohup ./generate_summary.sh chr18 > results/chr18.log &
# nohup ./generate_summary.sh chr19 > results/chr19.log &
# nohup ./generate_summary.sh chr20 > results/chr20.log &
# nohup ./generate_summary.sh chr21 > results/chr21.log &
# nohup ./generate_summary.sh chr22 > results/chr22.log &
#nohup ./generate_summary.sh chrX  > results/chrX.log &

echo -e "\nrunning"
