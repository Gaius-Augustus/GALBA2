#!/bin/bash
# Shared compute profile for all GALBA2 HPC test scenarios.
#
# This is the single source of truth for cluster/compute settings.
# Every scenario's run_test.sh sources this file first, then optionally
# sources scenario_overrides.sh for per-scenario tweaks.
#
# Override any variable ad-hoc:
#   PARTITION=highmem CORES=64 bash scenario_benchmark_athaliana_ep/run_test.sh

# Number of snakemake cores/jobs
CORES=${CORES:-48}

# SLURM partition(s), comma-separated
PARTITION=${PARTITION:-"batch,snowball,pinky"}

# Whether to use SLURM executor (set to false for local execution)
USE_SLURM=${USE_SLURM:-true}

# Per-job resources (exported as GALBA2_* env vars for config.ini override)
export GALBA2_CPUS_PER_TASK=${GALBA2_CPUS_PER_TASK:-48}
export GALBA2_MEM_OF_NODE=${GALBA2_MEM_OF_NODE:-120000}
export GALBA2_MAX_RUNTIME=${GALBA2_MAX_RUNTIME:-120}

# Default memory bound for SLURM --mem requests
DEFAULT_MEM_MB=${DEFAULT_MEM_MB:-120000}
