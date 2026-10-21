#!/bin/bash
# EpiSignaLR: methylation episignature classification for ONT, PacBio and
# Illumina 5-base data. Run with --help for usage.

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
source "$REPO_ROOT/lib/common.sh"
source "$REPO_ROOT/lib/steps.sh"

DEFAULT_MODELS="SETD1B,Kleefstra,Kabuki"

usage() {
    cat <<EOF
Usage: $(basename "$0") --platform <ont|pacbio|illumina5base> --sample <ID> --outdir <DIR>
                    (--input <BAM|gVCF> | --pileup <PILEUP>) [options]

Required:
  -p, --platform   Sequencing platform: ont, pacbio, illumina5base
  -s, --sample     Sample ID (used as output prefix and subdirectory name)
  -o, --outdir     Output directory; results are written to <outdir>/<sample>/
  -i, --input      Input file: modBAM (ont), aligned BAM (pacbio), or
                   methylation gVCF (illumina5base)
      --pileup     Use an existing pileup instead of --input and skip the
                   pileup step (modkit bedMethyl, pb-CpG-tools combined BED,
                   or gvcf_cpg_to_bed.py output, matching --platform)

Options:
  -r, --ref        hg38 reference FASTA (required for ont/pacbio with --input)
  -t, --threads    Threads for the pileup step (default: 4)
  -m, --models     Comma-separated LRS-derived models to run
                   (default: $DEFAULT_MODELS)
  -h, --help       Show this help and exit
EOF
}

PLATFORM=""
SAMPLE=""
OUT_ROOT=""
INPUT=""
PILEUP=""
REF=""
THREADS=4
MODELS_ARG="$DEFAULT_MODELS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--platform) PLATFORM="${2:-}"; shift 2 ;;
        -s|--sample)   SAMPLE="${2:-}"; shift 2 ;;
        -o|--outdir)   OUT_ROOT="${2:-}"; shift 2 ;;
        -i|--input)    INPUT="${2:-}"; shift 2 ;;
        --pileup)      PILEUP="${2:-}"; shift 2 ;;
        -r|--ref)      REF="${2:-}"; shift 2 ;;
        -t|--threads)  THREADS="${2:-}"; shift 2 ;;
        -m|--models)   MODELS_ARG="${2:-}"; shift 2 ;;
        -h|--help)     usage; exit 0 ;;
        *) usage >&2; die "Unknown argument: $1" ;;
    esac
done

# ---- Validate arguments -----------------------------------------------------

PLATFORM=$(echo "$PLATFORM" | tr '[:upper:]' '[:lower:]')
case "$PLATFORM" in
    ont|nanopore)                  PLATFORM=ont ;;
    pacbio|pb)                     PLATFORM=pacbio ;;
    illumina5base|ill5base|illumina) PLATFORM=illumina5base ;;
    "") usage >&2; die "--platform is required" ;;
    *)  die "Unsupported platform '$PLATFORM' (choose ont, pacbio, illumina5base)" ;;
esac

[[ -n "$SAMPLE" ]]   || { usage >&2; die "--sample is required"; }
[[ -n "$OUT_ROOT" ]] || { usage >&2; die "--outdir is required"; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer (got '$THREADS')"

if [[ -n "$INPUT" && -n "$PILEUP" ]]; then
    die "Specify only one of --input or --pileup"
elif [[ -z "$INPUT" && -z "$PILEUP" ]]; then
    usage >&2; die "One of --input or --pileup is required"
fi

source "$REPO_ROOT/lib/platforms/${PLATFORM}.sh"

if [[ -n "$INPUT" ]]; then
    [[ -f "$INPUT" ]] || die "Input file not found: $INPUT"
    INPUT=$(readlink -f "$INPUT")
    if [[ "$PLATFORM_REQUIRES_REF" -eq 1 ]]; then
        [[ -n "$REF" ]] || die "--ref is required for platform '$PLATFORM'"
        [[ -f "$REF" ]] || die "Reference not found: $REF"
        REF=$(readlink -f "$REF")
    fi
    require_cmd "${PLATFORM_DEPS[@]}"
else
    [[ -f "$PILEUP" ]] || die "Pileup file not found: $PILEUP"
    PILEUP=$(readlink -f "$PILEUP")
fi

IFS=',' read -r -a MODELS <<< "$MODELS_ARG"
for model in "${MODELS[@]}"; do
    [[ -d "$REPO_ROOT/data/${model}_enet_svm_model" ]] || \
        die "Unknown model '$model': $REPO_ROOT/data/${model}_enet_svm_model not found"
done

require_cmd python bedtools zcat

# ---- Set up output directories ------------------------------------------------

OUT_DIR="${OUT_ROOT%/}/${SAMPLE}"
PILEUP_DIR="$OUT_DIR/Pileups"
ARRAY_DIR="$OUT_DIR/Array_Derived_Episignature_Results"
LRS_DIR="$OUT_DIR/LRS_Derived_Episignature_Results"
ARRAY_CPG_BED="$REPO_ROOT/data/collapsed_coordinates_sorted.bed"

mkdir -p "$PILEUP_DIR" "$ARRAY_DIR" "$LRS_DIR"

log "EpiSignaLR: sample=$SAMPLE platform=$PLATFORM threads=$THREADS models=${MODELS[*]}"

# ---- Pileup -------------------------------------------------------------------

if [[ -n "$PILEUP" ]]; then
    log "Using provided pileup, skipping pileup step: $PILEUP"
else
    PILEUP=$(platform_pileup_path)
    if ! step_done "$PILEUP" "Pileup"; then
        log "Creating $PLATFORM pileup for $INPUT"
        platform_make_pileup "$INPUT"
        [[ -s "$PILEUP" ]] || die "Pileup step finished but expected output is missing: $PILEUP"
    fi
fi

# ---- Classification -------------------------------------------------------------

run_array_svm "$PILEUP"
run_lrs_models "$PILEUP"

log "EpiSignaLR finished for $SAMPLE. Results in $OUT_DIR"
