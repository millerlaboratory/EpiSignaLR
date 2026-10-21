#!/bin/bash
# Oxford Nanopore: modBAM -> modkit bedMethyl pileup.

PLATFORM_REQUIRES_REF=1
PLATFORM_DEPS=(modkit)
ARRAY_BED_TAG="modkit"

platform_pileup_path() {
    echo "${PILEUP_DIR}/${SAMPLE}.modkit.pileup_combined.bed.gz"
}

# platform_make_pileup <input_bam>
platform_make_pileup() {
    modkit pileup "$1" "$PILEUP_DIR" \
        --modified-bases 5mC 5hmC \
        --filter-threshold C:0.7 \
        --mod-threshold m:0.7 --mod-threshold h:0.7 \
        --cpg --combine-strands --phased --combine-mods \
        --ref "$REF" \
        --prefix "${SAMPLE}.modkit.pileup" --bgzf \
        -t "$THREADS" \
        --log-filepath "${PILEUP_DIR}/${SAMPLE}.log"
}

# stdin: decompressed pileup; stdout: 11-column bedMethyl-like rows, 'chr' stripped.
platform_to_array_bed() {
    awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | cut -f 1-11
}

# platform_region_average <samples_txt> <model_bed> <out_prefix>
platform_region_average() {
    python "$REPO_ROOT/scripts/region_average_methylation.py" \
        -s "$1" -b "$2" -o "$3" \
        -f 'csv' \
        --platform 'ont' \
        --min-coverage 2 \
        --weighted
}
