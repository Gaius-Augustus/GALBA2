#!/bin/bash
# Benchmark scenario: A. thaliana EP mode (proteins only, GALBA2)
# Data: /home/nas-hs/projs/braker4/data/Arabidopsis_thaliana/
#   - genome.fa (~121 MB, already softmasked)
#   - proteins.fa (~467 MB)
#   - ref_annot.gff3 (reference annotation for gffcompare)

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================================"
echo "Benchmark: A. thaliana — GALBA2 (protein evidence only)"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  - Genome:    /home/nas-hs/projs/braker4/data/Arabidopsis_thaliana/genome.fa"
echo "  - Proteins:  /home/nas-hs/projs/braker4/data/Arabidopsis_thaliana/proteins.fa (~467 MB)"
echo "  - Reference: ref_annot.gff3 (gffcompare evaluation)"
echo "  - Masking:   none (genome is pre-softmasked)"
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
snakemake \
    --snakefile "$PIPELINE_DIR/Snakefile" \
    --cores "$CORES" --jobs "$CORES" \
    $DRY_RUN_FLAG \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home,/home/nas-hs --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "Dry-run completed!" || echo "Benchmark completed!"
