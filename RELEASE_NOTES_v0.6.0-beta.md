# GALBA2 v0.6.0-beta — Release Notes

Release date: 2026-10-11
Previous release: [v0.5.0-beta](https://github.com/Gaius-Augustus/GALBA2/releases/tag/v0.5.0-beta)

Gene predictions are unchanged by this release: the toy scenario gives the
same gene set as v0.5.0-beta (0 differing lines in `galba.gtf` and `galba.aa`),
locally and through the SLURM executor on an HPC cluster.

## Node-local scratch

- Rules that write many or large temporary files (`run_augustus_hints`,
  `optimize_augustus`, `busco_genome`, `busco_proteins`, `run_compleasm`,
  `fantasia_annotate`) now work in a private directory on the compute node's
  local disk, copy only the kept files to `output/<sample>/`, and remove the
  scratch directory on exit, also on failure ([scripts/tmp_dir.sh](scripts/tmp_dir.sh),
  ported from BRAKER4 v0.6.0-beta). On network file systems with slow
  small-file writes this removes the dominant cost of these rules.
- The scratch root is `[paths] tmp_dir` in `config.ini` (or `GALBA2_TMP_DIR`);
  empty means each job's `$TMPDIR`, else `/tmp`. A root outside `/tmp` and
  `/home` is bound into the containers automatically. A job that finds less
  free space than it needs works in `output/<sample>/` as before and logs
  `WARNING: working in ...`; a job on scratch logs `scratch directory: ...`.
  Space requirements per rule are in the README ("Node-local scratch").
- `use_dev_shm = 1` now only selects `/dev/shm` as the scratch root of masking
  and the two AUGUSTUS rules.

## Bug fixes

- **compleasm failed on every run since the odb12.2 release of the BUSCO
  data** (compleasm issue #61): compleasm 0.2.7 crashed while parsing the
  current `file_versions.tsv` for placement files it only needs with
  `--autolineage`, and left a `placement_files.tmp` that made the next run
  stop with "another process is downloading". Two fixes: the
  `katharinahoff/galba2-tools` image now ships compleasm 0.2.10 (tag `v0.5.0`
  and `latest`), and [scripts/compleasm_wrapper.py](scripts/compleasm_wrapper.py)
  makes older images skip the placement download when a lineage is given.
  The wrapper also restores compleasm's offline fallback and skips proteins
  of 100,000 aa or more that hmmsearch cannot handle (BRAKER4 #91, #92, #105).
  A compleasm failure is now recorded as `COMPLEASM_FAILED: <reason>` in the
  summary and shown as a warning in the HTML report instead of passing as a
  result; BUSCO failures are shown the same way.
- **BUSCO proteome assessment** skips proteins of 100,000 aa or more, which
  made hmmsearch abort the whole run.
- **Gene lines wider than their transcripts** (BRAKER4 #62): the DIAMOND and
  transcript-ID filters re-span transcript and gene lines after dropping
  transcripts ([scripts/gtf_transcript_filter.py](scripts/gtf_transcript_filter.py));
  `normalize_gtf.py` splits gene IDs used on several loci and drops
  transcripts with an identical CDS chain (BRAKER4 #97, #71);
  [scripts/check_gtf_loci.py](scripts/check_gtf_loci.py) fails the pipeline
  on an inconsistent gene model after `normalize_cds` and on the final
  `galba.gtf`.
- **Infernal**: cmscan runs with `--oclan`, and `infernal_to_gff3.py` keeps
  the best hit of an overlapping set instead of discarding the whole locus
  (the old filter dropped every overlapping hit, rRNA clans included);
  column indices follow the Infernal 1.1.5 `--fmt 2` layout.
- **OMArk**: `[OMARK] ete_taxa_path` now works. `taxa.sqlite` is built once
  from `taxdump.tar.gz` under a lock and passed to OMArk with `-e`, so OMArk
  needs no internet on compute nodes.
- **gffcompare**: a gzipped `reference_gtf` is decompressed instead of
  silently yielding 0 CDS lines.
- **FANTASIA**: an empty `slurm_extra` made sbatch fail with "Unable to open
  file"; `--exclude` is only passed when `exclude_nodes` is set. A placeholder
  figure is written when no GO term passes the cutoff, so the rule no longer
  fails on an empty result. New rule `fantasia_gpu_check` probes the GPU at
  the start of a run and stops it within minutes if none is usable.
- **Snakemake**: `workflow.cores` is read defensively (it raises with
  `--executor slurm` and no `--cores`), with a warning when `--cores` caps
  `cpus_per_task`. Hard-coded 10 and 30 minute runtimes replaced.
- `compleasm` lineage handling: `-l NAME --odb odbNN` as given in
  `samples.csv` instead of forcing `_odb12` (BRAKER4 #73); `command -v`
  instead of `which` inside containers (issue #4); FANTASIA lookup bundle
  documented in `config.ini.example` (issue #3).
- `merge_ncrna_gff3.py` uses a bisection interval index; the merge was
  quadratic. Scripts referenced by rules but never committed
  (`merge_ncrna_gff3.py`, `rfam_family_types.tsv`, `generate_embeddings.py`)
  are now tracked.

## Output format

- `galba.gff3` follows NCBI/Ensembl conventions: coding transcripts are
  `mRNA` and genes carry `gene_biotype=protein_coding`; ncRNAs are written as
  gene -> RNA -> exon with `gene_biotype`, orphaned tRNA exons fixed and
  duplicate Infernal tRNAs dropped (Annotrieve compatibility).
- Fewer intermediate files are kept; BUSCO, compleasm and OMArk directories
  hold only their summaries after a run.

## Dependencies and infrastructure

- Host-side dependencies are managed with [pixi](https://pixi.sh)
  (`pixi.toml`, `pixi.lock`; `pixi install`, `pixi run test`).
  `requirements.txt` keeps the same pins for pip users; `pandas>=2.2.3,<3`
  as required by snakemake-executor-plugin-slurm 2.6.0.
- Container images: `katharinahoff/galba2-tools:v0.5.0` (= `latest`;
  compleasm 0.2.10, miniprot 0.14 first on PATH, minisplice),
  gffcompare 0.12.10, BUSCO 6.1.0. Build script
  [dockerfiles/build_push_galba2-tools.sh](dockerfiles/build_push_galba2-tools.sh).
- BUSCO citation updated to Tegenfeldt et al. 2025.
- New test scenario `test_scenarios/scenario_01_toy_slurm` (toy data through
  the SLURM executor with scratch, optimization, BUSCO, compleasm and
  FANTASIA) and 65 unit tests under `tests/`.

## Known limitations

- FANTASIA and OMArk were not re-run on a cluster for this release
  (the FANTASIA GPU job could not be scheduled in time); the code paths are
  shared with BRAKER4 v0.6.0-beta, where they are exercised.
