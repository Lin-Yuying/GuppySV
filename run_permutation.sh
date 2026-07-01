#!/usr/bin/env bash
set -euo pipefail

############################
# User-defined input files #
############################

INV="inv.max.bed"                 # inversion intervals (BED3)
SNP="${1}_filtered_positions.bed"       # SNP positions (BED3)
GENOME="female.guppy.genome.correct.sizes"         # chromosome sizes
BLACKLIST="blacklist.bed"
N=100000                       # number of permutations

########################
# Output files         #
########################

PERM_OUT="${1}_perm_results.txt"
LOG="${1}_permutation.log"

########################
# Start logging        #
########################

echo "Permutation test started at $(date)" > $LOG
echo "Inversion file: $INV" >> $LOG
echo "SNP file: $SNP" >> $LOG
echo "Genome file: $GENOME" >> $LOG
echo "Permutations: $N" >> $LOG
echo "-----------------------------" >> $LOG

########################
# Observed statistic   #
########################

OBS=$(bedtools intersect -u \
    -a "$INV" \
    -b "$SNP" \
    | wc -l)

echo "Observed SNPs in inversions: $OBS" | tee -a $LOG

########################
# Permutation test     #
########################

> "$PERM_OUT"

for i in $(seq 1 $N); do
    bedtools shuffle \
        -i "$SNP" \
        -g "$GENOME" \
        -incl allowed_regions.bed \
    | bedtools intersect -u \
        -a "$INV" \
        -b - \
    | wc -l >> "$PERM_OUT"

    # Progress update every 1000 iterations
    if (( i % 1000 == 0 )); then
        echo "Completed $i permutations" >> $LOG
    fi
done

########################
# Empirical p-value    #
########################

PVAL=$(awk -v obs="$OBS" '
    $1 >= obs {count++}
    END {print (count+1)/(NR+1)}
' "$PERM_OUT")

NULL_MEAN=$(awk '{sum+=$1} END {print sum/NR}' "$PERM_OUT")

########################
# Final report         #
########################

echo "-----------------------------" >> $LOG
echo "Null mean overlap: $NULL_MEAN" | tee -a $LOG
echo "Empirical p-value: $PVAL" | tee -a $LOG
echo "Permutation test finished at $(date)" >> $LOG

########################
# Print summary        #
########################

echo
echo "===== PERMUTATION TEST SUMMARY ====="
echo "Observed SNPs in inversions : $OBS"
echo "Null mean (10,000 perms)    : $NULL_MEAN"
echo "Empirical p-value           : $PVAL"
echo "==================================="

