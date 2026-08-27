CHR=$1
START=$2
END=$3
pixy --stats dxy pi --vcf all.$CHR.vcf.gz --populations ${CHR}_inversion.popfile.txt --output_folder ./pixy_output --output_prefix $CHR.inv.10k  --bypass_invariant_check --window_size 10000 --interval_start ${START} --interval_end ${END}&
