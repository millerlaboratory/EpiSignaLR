#!/bin/bash
# Platform-independent classification steps. Sourced by EpiSignaLR.sh.
# Expects SAMPLE, ARRAY_DIR, LRS_DIR, PILEUP_DIR, REPO_ROOT, MODELS and the
# platform module functions to be defined.

# run_array_svm <pileup>
# Intersect the pileup with the array-derived episignature CpGs and run the SVMs.
run_array_svm() {
    local pileup=$1
    local array_bed="${ARRAY_DIR}/${SAMPLE}.${ARRAY_BED_TAG}.episignature.bed"

    if ! step_done "$array_bed" "Array episignature BED"; then
        log "Intersecting pileup with array episignature CpGs for $SAMPLE..."
        # Write to a temp file so a failed run is not mistaken for a finished one.
        zcat "$pileup" | platform_to_array_bed | \
            bedtools intersect -a "$ARRAY_CPG_BED" -b - -loj -wa -wb > "${array_bed}.tmp"
        mv "${array_bed}.tmp" "$array_bed"
    fi

    if ! step_done "${ARRAY_DIR}/${SAMPLE}.results_scores.csv" "Array-derived SVM results"; then
        log "Running Array Derived SVMs for $SAMPLE..."
        python "$REPO_ROOT/scripts/run_SVM.py" \
            --bed "$array_bed" \
            --sample "$SAMPLE" \
            --out_prefix "${ARRAY_DIR}/${SAMPLE}.results"
    fi

    log "Finished running Array Derived SVMs for $SAMPLE"
}

# run_lrs_models <pileup>
# Region-average methylation + ElasticNet/SVM for each LRS-derived model in MODELS.
run_lrs_models() {
    local pileup=$1
    local samples_txt="${PILEUP_DIR}/${SAMPLE}.samples.txt"
    local model model_dir region_csv

    # Second column sets the sample name used in the region matrix and results.
    printf '%s\t%s\n' "$pileup" "$SAMPLE" > "$samples_txt"

    for model in "${MODELS[@]}"; do
        model_dir="$REPO_ROOT/data/${model}_enet_svm_model"
        region_csv="${LRS_DIR}/${SAMPLE}.${model}_enet_svm_full_panel_DMPs.csv"

        log "Running LRS Derived ElasticNet+SVM ($model) for $SAMPLE..."

        if ! step_done "$region_csv" "$model region average methylation"; then
            platform_region_average \
                "$samples_txt" \
                "${model_dir}/${model}_enet_svm_full_panel_DMPs.bed" \
                "${LRS_DIR}/${SAMPLE}.${model}_enet_svm_full_panel_DMPs"
        fi

        if ! step_done "${LRS_DIR}/${SAMPLE}.${model}.svm_results_top_result.txt" "$model SVM results"; then
            python "$REPO_ROOT/scripts/run_enet_svm.py" \
                --region-matrix "$region_csv" \
                --model-dir "${model_dir}/" \
                --out-prefix "${LRS_DIR}/${SAMPLE}.${model}.svm_results"
        fi
    done

    log "Finished running LRS Derived ElasticNet+SVMs for $SAMPLE"
}
