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

if [ -f "${OUT_DIR}/${prefix}.modkit.pileup_combined.bed.gz" ]; then
    echo "Pileup already exists, skipping pileup step"
else
    echo "Creating pileup for $BAM"

    modkit pileup ${BAM} "${OUT_DIR}" \
    --modified-bases 5mC 5hmC \
    --filter-threshold C:0.7 \
    --mod-threshold m:0.7 --mod-threshold h:0.7 \
    --cpg --combine-strands --phased --combine-mods \
    --ref $REF \
    --prefix ${prefix}.modkit.pileup --bgzf \
    -t 30 \
    --log-filepath ${OUT_DIR}/${prefix}.log
fi

zcat "${OUT_DIR}/${prefix}.modkit.pileup_combined.bed.gz" | awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | cut -f 1-11 | \
    bedtools intersect -a "$BED" -b - -loj -wa -wb > "${OUT_DIR}/${prefix}.modkit.episignature.bed"

echo "Finished processing $BAM"


python "$REPO_ROOT/scripts/run_SVM.py" \
    --bed "${OUT_DIR}/${prefix}.modkit.episignature.bed" \
    --out_prefix "${OUT_DIR}/${prefix}.results"

