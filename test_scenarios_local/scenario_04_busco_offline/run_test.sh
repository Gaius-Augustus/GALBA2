#!/bin/bash
# Scenario 04: BUSCO offline mode — busco_download_path points at a
# pre-downloaded lineage (with dataset.cfg), so BUSCO should run --offline
# and compleasm should reuse the same lineages/ directory.
set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

source "$TIER_DIR/compute_profile.sh"
export GALBA2_CONFIG="$SCENARIO_DIR/config.ini"

cd "$SCENARIO_DIR"
export SINGULARITYENV_PREPEND_PATH=/opt/conda/bin
snakemake \
    --snakefile "$PIPELINE_DIR/Snakefile" \
    --cores "$CORES" --jobs "$CORES" \
    --printshellcmds --rerun-incomplete --nolock \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS
