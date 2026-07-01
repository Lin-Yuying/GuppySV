chr=$1

prefix="all"

vcftools --gzvcf ../VCF/${prefix}.LG${chr}.vcf.gz --maf 0.05 --thin 100  -c --geno-r2 --max-missing-count 0 | perl /Users/evolutioneco/Project/3.SV/4.localPCA/emerald2windowldcounts.pl | gzip >${prefix}.thin100.maf5.LG1.windows.ld.gz

