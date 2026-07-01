prefix=$1

# convert vcf to other formats using plink2
plink2 --vcf ${prefix}.filter.recode.vcf.gz --make-bed --out ${prefix} --allow-extra-chr

# PCA 
plink2 --bfile ${prefix} --pca --out ${prefix} --allow-extra-chr

# admixture
for k in $(seq 1 8)
do
    ./admixture --cv ${prefix}.bed ${k} | tee log${k}.out
done
