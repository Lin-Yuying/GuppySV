#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <CHROM> <START> <END>"
    exit 1
fi

CHROM=$1
START=$2
END=$3

VCF="../all.${CHROM}.vcf.gz"

# ===== females =====
bcftools view \
  -S females.txt \
  -r ${CHROM}:${START}-${END} \
  -v snps \
  $VCF \
| bcftools query -f '%CHROM\t%POS[\t%DP]\n' \
| awk '{
    sum=0;
    for(i=3;i<=NF;i++) sum+=$i;
    mean=sum/(NF-2);
    print $1":"$2"\t"mean
}' > ${CHROM}_${START}_${END}_female.depth

# ===== males =====
bcftools view \
  -S males.txt \
  -r ${CHROM}:${START}-${END} \
  -v snps \
  $VCF \
| bcftools query -f '%CHROM\t%POS[\t%DP]\n' \
| awk '{
    sum=0;
    for(i=3;i<=NF;i++) sum+=$i;
    mean=sum/(NF-2);
    print $1":"$2"\t"mean
}'  > ${CHROM}_${START}_${END}_male.depth

# ===== merge =====
join -1 1 -2 1 ${CHROM}_${START}_${END}_female.depth ${CHROM}_${START}_${END}_male.depth \
| awk -F'[:\t]' '{
    print $1"\t"$2"\t"$3"\t"$4
}' > ${CHROM}_${START}_${END}_MF.depth 

ln -s ${CHROM}_${START}_${END}_MF.depth ${CHROM}.MF.depth
