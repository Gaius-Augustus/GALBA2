#!/bin/bash
# Per-scenario overrides for the toy data run through the SLURM executor.
#
# Smoke test of the whole DAG on the cluster: every rule converted to
# node-local scratch (run_augustus_hints, optimize_augustus, busco_genome,
# busco_proteins, run_compleasm, fantasia_annotate) runs on real compute
# nodes with [paths] tmp_dir empty, i.e. the job's $TMPDIR (/tmp on brain).
# The toy genome (~1 Mb) keeps every job short; optimize_augustus.pl is ON
# so that its scratch path is exercised (2-fold with < 400 training genes).

# Small jobs: 16 CPUs per job instead of the 48 of the benchmark profile.
export CORES=16
export GALBA2_CPUS_PER_TASK=16
export GALBA2_SKIP_OPTIMIZE_AUGUSTUS=0
export GALBA2_SKIP_BUSCO=0
export GALBA2_NO_CLEANUP=1
export GALBA2_USE_DEV_SHM=0
export GALBA2_MAX_RUNTIME=240
export DEFAULT_MEM_MB=60000
export GALBA2_MEM_OF_NODE=60000

# brain HPC: Snakemake is in hoffk83's conda env, not on the default PATH.
export SNAKEMAKE_BIN=/home/hoffk83/miniconda3/bin/snakemake

# brain HPC: the FANTASIA image and lookup bundle live under /projects.
export SINGULARITY_EXTRA_BIND=",/projects"

# brain HPC: singularity is a module (module avail lists singularity/3.10.0
# and singularity/3.11.3); load it so snakemake and the SLURM jobs, which
# inherit this PATH, can find it.
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity/3.11.3 2>/dev/null || true
