#!/bin/bash
# Per-scenario overrides for toy test.
# Skip optimization for fast testing (~5 min instead of ~30 min).

export GALBA2_SKIP_OPTIMIZE_AUGUSTUS=1
export GALBA2_SKIP_BUSCO=1
export GALBA2_NO_CLEANUP=1
