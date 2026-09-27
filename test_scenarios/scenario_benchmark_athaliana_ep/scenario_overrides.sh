#!/bin/bash
# Per-scenario overrides for A. thaliana GALBA2 benchmark.
#
# Full A. thaliana genome with Viridiplantae proteins.
# Runs optimize_augustus.pl (the slow but proper way).

export GALBA2_USE_DEV_SHM=1
export GALBA2_NO_CLEANUP=1
export GALBA2_MAX_RUNTIME=4320
export DEFAULT_MEM_MB=250000

# On: run optimize_augustus.pl (matches original galba.pl behavior).
export GALBA2_SKIP_OPTIMIZE_AUGUSTUS=0

# brain HPC: Snakemake is in hoffk83's conda env, not on the default PATH.
export SNAKEMAKE_BIN=/home/hoffk83/miniconda3/bin/snakemake

# brain HPC: genome and reference annotation are under /projects — bind-mount it.
export SINGULARITY_EXTRA_BIND=",/projects"

# brain HPC: singularity is a module; load it so snakemake can find it.
# The || true prevents set -e from aborting if the module system is unavailable.
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity 2>/dev/null || true
