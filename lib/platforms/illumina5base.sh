#!/bin/bash
# Illumina 5-base: DRAGEN methylation gVCF -> per-CpG BED.

PLATFORM_REQUIRES_REF=0
PLATFORM_DEPS=()
ARRAY_BED_TAG="array"

platform_pileup_path() {
    echo "${PILEUP_DIR}/${SAMPLE}_ill5base_cpg.bed.gz"
}

# platform_make_pileup <input_gvcf>
platform_make_pileup() {
    python "$REPO_ROOT/scripts/gvcf_cpg_to_bed.py" \
        --gvcf "$1" \
        --output "$(platform_pileup_path)" \
        --no-header
}

# stdin: decompressed pileup; stdout: 11-column bedMethyl-like rows, 'chr' stripped.
platform_to_array_bed() {
    awk 'BEGIN {FS=OFS="\t"} {$1 = substr($1, 4)} 1' | \
        awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"m",$5,".",$2,$3,".",$5,$4*100}'
}

# platform_region_average <samples_txt> <model_bed> <out_prefix>
platform_region_average() {
    python "$REPO_ROOT/scripts/region_average_methylation_cpgbed.py" \
        -s "$1" -b "$2" -o "$3" \
        -f 'csv' \
        --weighted
}
