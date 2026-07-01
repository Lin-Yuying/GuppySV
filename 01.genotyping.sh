bcftools=~/bin/bcftools-1.16/bcftools
ref=~/3.SV/1.alignment_male/guppy_female.fa
region_name=$1 #batch calling, I seperate it into 4 regions
bam_file=$2 #please change when use
${bcftools} mpileup -Ou -q 20 -Q 20 -R region${region_name}.guppy.female.txt --skip-indels -a FORMAT/AD,FORMAT/DP -f $ref -b ${bam_file}.pop.txt | ${bcftools} call -mv -Oz -a GQ,GP -o ${bam_file}.region${region_name}.vcf.gz

# SNP filtering
input=$1
output=$2
vcftools --gzvcf ${input} \
--maf 0.01 \
--min-alleles 2 \
--max-alleles 2 \
--max-missing 0.9 \
--min-meanDP 10 \
--max-meanDP 80 \
--max-missing 0.8 \
--recode \
--recode-INFO-all \
--minGQ 20 \
--minQ 200 \
--out ${output}

bgzip ${output}.recode.vcf
