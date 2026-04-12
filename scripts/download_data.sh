#!/bin/bash
# Download shared databases for the GALBA2 pipeline.
#
# Usage:
#   bash scripts/download_data.sh                        # auto-detect location
#   bash scripts/download_data.sh --shared-data /path    # custom location
#   bash scripts/download_data.sh --omark                # also download LUCA.h5
#   bash scripts/download_data.sh --all                  # download everything
#
# Databases:
#   - compleasm lineage libraries (downloaded on-demand by compleasm, but
#     pre-downloading avoids repeated fetches in multi-sample runs)
#   - BUSCO lineage datasets (same rationale)
#   - OMAmer LUCA.h5 (~15 GB, only with --omark or --all)
#
# Shared data location resolution:
#   1. If --shared-data <path> is given, use that
#   2. If shared_data/ exists in this repo, use it
#   3. If a sibling BRAKER4 or BRAKER-as-snakemake repo has shared_data/, reuse it
#   4. Otherwise, create shared_data/ in this repo
#
# The resolved path is printed at the end so you can set it in config.ini:
#   [paths]
#   busco_download_path = <resolved>/busco_downloads
#   compleasm_download_path = <resolved>/compleasm_downloads

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Defaults
SHARED_DATA=""
DOWNLOAD_OMARK=0
DOWNLOAD_ALL=0

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --shared-data)
            shift_next=1
            ;;
        --omark)
            DOWNLOAD_OMARK=1
            ;;
        --all)
            DOWNLOAD_ALL=1
            DOWNLOAD_OMARK=1
            ;;
        -h|--help)
            sed -n '2,28p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            if [ "${shift_next:-0}" = "1" ]; then
                SHARED_DATA="$arg"
                shift_next=0
            else
                echo "Unknown argument: $arg (use -h for help)" >&2
                exit 1
            fi
            ;;
    esac
done

# ============================================================================
# Resolve shared_data location
# ============================================================================
resolve_shared_data() {
    # 1. Explicit --shared-data
    if [ -n "$SHARED_DATA" ]; then
        echo "$SHARED_DATA"
        return
    fi

    # 2. This repo's own shared_data
    if [ -d "$REPO_DIR/shared_data" ]; then
        # Check if it has actual content (not just empty dirs)
        if [ -d "$REPO_DIR/shared_data/compleasm_downloads" ] || \
           [ -d "$REPO_DIR/shared_data/busco_downloads" ]; then
            echo "$REPO_DIR/shared_data"
            return
        fi
    fi

    # 3. Sibling repos (one level up)
    PARENT_DIR="$(dirname "$REPO_DIR")"
    for sibling in BRAKER4 BRAKER-as-snakemake; do
        SIBLING_SHARED="$PARENT_DIR/$sibling/shared_data"
        if [ -d "$SIBLING_SHARED" ]; then
            if [ -d "$SIBLING_SHARED/compleasm_downloads" ] || \
               [ -d "$SIBLING_SHARED/busco_downloads" ]; then
                echo "  Found existing shared_data in sibling $sibling repo" >&2
                echo "$SIBLING_SHARED"
                return
            fi
        fi
    done

    # 4. Fallback: create in this repo
    echo "$REPO_DIR/shared_data"
}

SHARED_DATA_DIR="$(resolve_shared_data)"
SHARED_DATA_DIR="$(cd "$REPO_DIR" && mkdir -p "$SHARED_DATA_DIR" && cd "$SHARED_DATA_DIR" && pwd)"

echo "============================================================"
echo "GALBA2 Shared Data Downloader"
echo "============================================================"
echo ""
echo "Shared data directory: $SHARED_DATA_DIR"
echo ""

# ============================================================================
# compleasm lineage libraries
# ============================================================================
COMPLEASM_DIR="$SHARED_DATA_DIR/compleasm_downloads"
mkdir -p "$COMPLEASM_DIR"

if [ -d "$COMPLEASM_DIR" ] && [ "$(ls -A "$COMPLEASM_DIR" 2>/dev/null)" ]; then
    N_LINEAGES=$(ls -d "$COMPLEASM_DIR"/mb_downloads/*/lineages/* 2>/dev/null | wc -l || echo 0)
    echo "compleasm lineage libraries: $N_LINEAGES lineage(s) already cached in $COMPLEASM_DIR"
    echo "  (compleasm downloads additional lineages on-demand during pipeline runs)"
else
    echo "compleasm lineage libraries: $COMPLEASM_DIR (empty — will be populated on first run)"
fi
echo ""

# ============================================================================
# BUSCO lineage datasets
# ============================================================================
BUSCO_DIR="$SHARED_DATA_DIR/busco_downloads"
mkdir -p "$BUSCO_DIR"

if [ -d "$BUSCO_DIR" ] && [ "$(ls -A "$BUSCO_DIR" 2>/dev/null)" ]; then
    N_BUSCO=$(ls -d "$BUSCO_DIR"/lineages/* 2>/dev/null | wc -l || echo 0)
    echo "BUSCO lineage datasets: $N_BUSCO lineage(s) already cached in $BUSCO_DIR"
    echo "  (BUSCO downloads additional lineages on-demand during pipeline runs)"
else
    echo "BUSCO lineage datasets: $BUSCO_DIR (empty — will be populated on first run)"
fi
echo ""

# ============================================================================
# OMAmer LUCA.h5 database (optional, ~15 GB; only for run_omark=1)
# ============================================================================
if [ "$DOWNLOAD_OMARK" = "1" ]; then
    LUCA_DB="$SHARED_DATA_DIR/LUCA.h5"

    # Also check sibling repos' test_data for existing LUCA.h5
    if [ ! -f "$LUCA_DB" ]; then
        PARENT_DIR="$(dirname "$REPO_DIR")"
        for sibling in BRAKER4 BRAKER-as-snakemake; do
            SIBLING_LUCA="$PARENT_DIR/$sibling/test_data/LUCA.h5"
            if [ -f "$SIBLING_LUCA" ]; then
                echo "Found LUCA.h5 in sibling $sibling repo: $SIBLING_LUCA"
                echo "  Creating symlink instead of re-downloading..."
                ln -sf "$SIBLING_LUCA" "$LUCA_DB"
                break
            fi
        done
    fi

    if [ ! -f "$LUCA_DB" ] && [ ! -L "$LUCA_DB" ]; then
        echo "Downloading OMAmer LUCA.h5 database (~15 GB)..."
        echo "  URL: https://omabrowser.org/All/LUCA.h5"
        wget -q --show-progress "https://omabrowser.org/All/LUCA.h5" -O "$LUCA_DB"
        echo "  Downloaded: $LUCA_DB ($(du -h "$LUCA_DB" | cut -f1))"
    else
        LUCA_SIZE=$(du -h "$LUCA_DB" | cut -f1)
        echo "OMAmer LUCA.h5 database: $LUCA_DB ($LUCA_SIZE)"
    fi
    echo ""
else
    echo "OMAmer LUCA.h5: skipped (use --omark or --all to download)"
    echo ""
fi

# ============================================================================
# Summary and config.ini guidance
# ============================================================================
echo "============================================================"
echo "Shared data ready at: $SHARED_DATA_DIR"
echo ""
echo "To use this location, add to your config.ini [paths] section:"
echo ""
echo "  [paths]"
echo "  busco_download_path = $BUSCO_DIR"
echo "  compleasm_download_path = $COMPLEASM_DIR"
if [ "$DOWNLOAD_OMARK" = "1" ] && [ -f "$SHARED_DATA_DIR/LUCA.h5" -o -L "$SHARED_DATA_DIR/LUCA.h5" ]; then
    LUCA_RESOLVED="$(readlink -f "$SHARED_DATA_DIR/LUCA.h5")"
    echo ""
    echo "For OMArk, also set omamer_db in config.ini [PARAMS]:"
    echo "  omamer_db = $LUCA_RESOLVED"
fi
echo ""
echo "Or if this is the default location (shared_data/ in this repo"
echo "or a sibling), no config change is needed — the pipeline"
echo "auto-discovers it."
echo "============================================================"
