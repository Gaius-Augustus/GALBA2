#!/bin/bash
# Shared compute profile for local (no-SLURM) test scenarios.
#
# This file is the SINGLE SOURCE OF TRUTH for the cores/memory configuration
# used by every scenario in test_scenarios_local/. Edit this file once to
# scale the local tests up or down.
#
# Every variable uses the ${VAR:-default} idiom so that an external
# environment variable still wins:
#   CORES=4 bash scenario_01_ep/run_test.sh

# How many cores snakemake itself uses (snakemake --cores / --jobs).
export CORES="${CORES:-8}"

# SLURM partition list — only consulted if USE_SLURM=true.
export PARTITION="${PARTITION:-batch}"

# Local mode: do NOT submit jobs to SLURM.
export USE_SLURM="${USE_SLURM:-false}"

# Per-job resource requests (still consulted by the Snakefile in local mode).
export GALBA2_CPUS_PER_TASK="${GALBA2_CPUS_PER_TASK:-8}"
export GALBA2_MEM_OF_NODE="${GALBA2_MEM_OF_NODE:-32000}"
export GALBA2_MAX_RUNTIME="${GALBA2_MAX_RUNTIME:-1440}"

# Maximum mem_mb that snakemake will request as a default-resource per job.
export DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-32000}"

# Build the snakemake --executor / --default-resources arguments.
if [ "$USE_SLURM" = "true" ]; then
    export EXECUTOR_ARGS="--executor slurm --default-resources slurm_partition=$PARTITION mem_mb=$DEFAULT_MEM_MB"
else
    export EXECUTOR_ARGS=""
fi
