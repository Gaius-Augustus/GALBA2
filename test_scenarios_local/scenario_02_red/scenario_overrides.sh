#!/bin/bash
# Per-scenario overrides for Red masking test.
# Leaves genome_masked empty so the masking rule runs; selects Red (fast).

export GALBA2_SKIP_OPTIMIZE_AUGUSTUS=1
export GALBA2_SKIP_BUSCO=1
export GALBA2_NO_CLEANUP=1
export GALBA2_MASKING_TOOL=red
