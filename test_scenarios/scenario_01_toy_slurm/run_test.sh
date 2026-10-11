#!/bin/bash
# Smoke test scenario: toy data (test_data/) through the SLURM executor

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================================"
echo "Smoke test: toy data through the SLURM executor (node-local scratch)"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  - Genome:    test_data/genome.fa (~1 MB)"
echo "  - Proteins:  test_data/proteins.fa (~119 KB)"
echo "  - Scratch:   [paths] tmp_dir empty -> each job uses its \$TMPDIR"
echo "  - Masking:   none (genome_masked set)"
echo ""

source "$TIER_DIR/compute_profile.sh"
[ -f "$SCENARIO_DIR/scenario_overrides.sh" ] && source "$SCENARIO_DIR/scenario_overrides.sh"

if [ "$USE_SLURM" = "true" ]; then
    EXECUTOR_ARGS="--executor slurm --default-resources slurm_partition=$PARTITION mem_mb=$DEFAULT_MEM_MB"
else
    EXECUTOR_ARGS=""
fi

export GALBA2_CONFIG="$TIER_DIR/config.ini"

DRY_RUN=${DRY_RUN:-false}
if [ "$DRY_RUN" = "true" ]; then
    DRY_RUN_FLAG="-n"
    echo "Running in DRY-RUN mode"
else
    DRY_RUN_FLAG=""
    echo "Running FULL EXECUTION"
fi

cd "$SCENARIO_DIR"

export SINGULARITYENV_PREPEND_PATH=/opt/conda/bin
# Allow scenario_overrides.sh to inject additional Singularity bind paths
# (e.g. SINGULARITY_EXTRA_BIND=",/projects" when input data lives under /projects).
SINGULARITY_EXTRA_BIND="${SINGULARITY_EXTRA_BIND:-""}"
$SNAKEMAKE_BIN \
    --snakefile "$PIPELINE_DIR/Snakefile" \
    --cores "$CORES" --jobs "$CORES" \
    $DRY_RUN_FLAG \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home${SINGULARITY_EXTRA_BIND} --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "Dry-run completed!" || echo "Smoke test completed!"
