#!/bin/bash

BAM=$1
prefix=$2
OUT_DIR=$3/${prefix}
mkdir -p $OUT_DIR

BED=.../data/collapsed_coordinates_sorted.bed

cd $OUT_DIR

module load modkit/0.4.4-rc1
module load samtools/1.22
module load bedtools/2.31.1

echo "Processing $BAM"

if [ -f "${OUT_DIR}/${prefix}.modkit.pileup.bed.gz" ]; then
    echo "Pileup already exists, skipping pileup step"
else
    echo "Creating pileup for $BAM"
    modkit pileup $BAM \
        "${OUT_DIR}/${prefix}.modkit.pileup.bed" \
        -t 20 \
        --ref #ADD REFERENCE FASTA FILE PATH \
        --preset traditional
    bgzip "${OUT_DIR}/${prefix}.modkit.pileup.bed"
fi

zcat "${OUT_DIR}/${prefix}.modkit.pileup.bed.gz" | awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | cut -f 1-11 | \
    bedtools intersect -a "$BED" -b - -loj -wa -wb > "${OUT_DIR}/${prefix}.modkit.episignature.bed"

echo "Finished processing $BAM"

source activate charge-kabuki-classifier

python ../scripts/run_SVM.py \
    --bed "${OUT_DIR}/${prefix}.modkit.episignature.bed" \
    --out_prefix "${OUT_DIR}/${prefix}.results"

