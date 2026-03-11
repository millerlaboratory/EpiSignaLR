#!/bin/bash

BAM=$1
prefix=$2
OUT_DIR=$3/${prefix}
REF=$4

mkdir -p $OUT_DIR

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

BED="$REPO_ROOT/data/collapsed_coordinates_sorted.bed"

cd $OUT_DIR

echo "Processing $BAM"

if [ -f "${OUT_DIR}/${prefix}.modkit.pileup.bed.gz" ]; then
    echo "Pileup already exists, skipping pileup step"
else
    echo "Creating pileup for $BAM"
    modkit pileup $BAM \
        "${OUT_DIR}/${prefix}.modkit.pileup.bed" \
        -t 20 \
        --ref "${REF}" \
        --preset traditional
    bgzip "${OUT_DIR}/${prefix}.modkit.pileup.bed"
fi

zcat "${OUT_DIR}/${prefix}.modkit.pileup.bed.gz" | awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | cut -f 1-11 | \
    bedtools intersect -a "$BED" -b - -loj -wa -wb > "${OUT_DIR}/${prefix}.modkit.episignature.bed"

echo "Finished processing $BAM"


python "$REPO_ROOT/scripts/run_SVM.py" \
    --bed "${OUT_DIR}/${prefix}.modkit.episignature.bed" \
    --out_prefix "${OUT_DIR}/${prefix}.results"

