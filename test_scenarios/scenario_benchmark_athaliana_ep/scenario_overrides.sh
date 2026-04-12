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
