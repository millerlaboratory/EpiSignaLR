#!/bin/bash
# Shared helpers for EpiSignaLR. Sourced by EpiSignaLR.sh.

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

require_cmd() {
    local cmd
    for cmd in "$@"; do
        command -v "$cmd" >/dev/null 2>&1 || die "Required command '$cmd' not found in PATH"
    done
}

# step_done <output_file> <description>
# Returns 0 (and logs a skip message) if the output already exists and is non-empty.
step_done() {
    if [[ -s "$1" ]]; then
        log "$2 already exists, skipping: $1"
        return 0
    fi
    return 1
}
