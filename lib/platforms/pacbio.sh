#!/bin/bash
# PacBio: aligned BAM -> pb-CpG-tools pileup.

PLATFORM_REQUIRES_REF=1
PLATFORM_DEPS=(aligned_bam_to_cpg_scores)
ARRAY_BED_TAG="pbcpg"

platform_pileup_path() {
    echo "${PILEUP_DIR}/${SAMPLE}_pbcpg.combined.bed.gz"
}

# platform_make_pileup <input_bam>
platform_make_pileup() {
    aligned_bam_to_cpg_scores \
        --bam "$1" \
        --ref "$REF" \
        --pileup-mode count \
        --modsites-mode reference \
        --output-prefix "${PILEUP_DIR}/${SAMPLE}_pbcpg" \
        --threads "$THREADS"
}

# stdin: decompressed pileup; stdout: 11-column bedMethyl-like rows, 'chr' stripped.
platform_to_array_bed() {
    sed '/^#/d' | \
        awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | \
        cut -f 1,2,3,4,6 | \
        awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"m",$5,".",$2,$3,".",$5,$4}'
}

# platform_region_average <samples_txt> <model_bed> <out_prefix>
platform_region_average() {
    python "$REPO_ROOT/scripts/region_average_methylation.py" \
        -s "$1" -b "$2" -o "$3" \
        -f 'csv' \
        --platform pacbio \
        --min-coverage 2 \
        --weighted
}
