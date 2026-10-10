# GALBA2 v0.4.0-beta — Release Notes

Release date: 2026-04-15
Previous release: [v0.3.0-beta](https://github.com/Gaius-Augustus/GALBA2/releases/tag/v0.3.0-beta)

## Bug fixes (accuracy-critical)

- **AUGUSTUS now runs with the correct hints against the masked genome.** The `splitMfasta` post-split rename loop in [rules/augustus_predict/run_augustus_hints.smk](rules/augustus_predict/run_augustus_hints.smk) had a hardcoded `genome.split.*` glob that silently failed when the input was `genome_masked.fa`, leaving chunks named by index rather than chromosome. AUGUSTUS could not match hints against sequence names, fell back to ab initio prediction, and locus Sn on *A. thaliana* dropped from ~70 to ~58. Fixed by widening the glob to `*.split.*`. (`026a390`)
- **Serialise OMArk via flock.** Parallel multi-sample OMArk runs were racing on the CWD-local `.etetoolkit/taxa.sqlite`, producing `sqlite3.OperationalError: database is locked`. The `omark` invocation is now wrapped in `flock` on a shared lockfile. (`166871e`)
- **OMArk wiring.** `omamer_db` now reads from the `[OMARK]` section of `config.ini` (matching BRAKER4 convention) and from the `GALBA2_OMAMER_DB` environment variable, with a `[PARAMS] omamer_db` back-compat fallback. Hardened `run_omark.smk` against `None`/`False` sentinels. (`9de7d33`)

- **Gene lines in `galba.gtf` now match their transcripts (BRAKER4 issue #62).** `filter_gtf_by_txid.py` and `filter_gtf_by_diamond_against_ref.py` (both v1.1.0) previously copied gene lines unchanged after dropping transcripts, so a gene line could be wider than its remaining transcripts (3 of 30549 genes in an *A. thaliana* run). Both scripts now use the shared module [scripts/gtf_transcript_filter.py](scripts/gtf_transcript_filter.py), which re-spans transcript lines to their features and gene lines to their remaining transcripts, and drops genes left without transcripts. Rule `sanity_check_augustus` additionally runs [scripts/check_gtf_loci.py](scripts/check_gtf_loci.py) on the final `galba.gtf` and fails on any gene/transcript span inconsistency; before, this check ran only on `galba.normalized.gtf` in rule `normalize_cds`.

## New features

- **Offline BUSCO / compleasm support** for firewalled clusters. BUSCO rules auto-detect `{download_path}/lineages/{lineage}/dataset.cfg` and switch to `--offline`. `compleasm_download_path` now defaults to `{busco_download_path}/lineages` so a single extracted lineage tarball serves both tools. A new README section explains the pre-download workflow. (`9cdecc8`)
- **`requirements.txt`** for host-side Python deps (`snakemake==8.18.2`, `snakemake-executor-plugin-slurm==2.6.0`, `pandas==2.0.2`). Enables GitHub dependency graph + Dependabot. (`751c6ad`)

## Housekeeping

- Removed the broken BUSCO-container plot path — BUSCO v6 no longer ships `generate_plot.py`/R/ggplot2. The HTML report's BUSCO+compleasm figure is generated in-house by [scripts/completeness_plot.py](scripts/completeness_plot.py) (matplotlib only, no external deps). (`ba322dc`)
- Trimmed README redundancy: merged the overlapping "What is different in GALBA2?" and "What is GALBA?" sections into one overview; dropped misleading "EP mode" naming (GALBA has only one prediction mode). Renamed `test_scenarios_local/scenario_01_ep/` → `scenario_01_toy/`. (`379c11c`)
- Benchmark scenarios now declare the pre-softmasked *A. thaliana* genome as pre-masked instead of leaving `genome_masked` empty — avoids unnecessary RepeatMasker runs now that masking is on the gene-prediction critical path. (`f3910c2`)
- README benchmark section rewritten with `galba.pl` vs GALBA2 as the primary comparison (the section title's long-standing promise), F1 scores added, bolding consistent with the "best per column" legend. Key finding: GALBA2 and `galba.pl` are accuracy-equivalent on *A. thaliana* (locus F1 66.7 vs 66.5, base F1 89.1 vs 89.1) — GALBA2's value is the Snakemake execution model (multi-sample input, automatic resume, SLURM-native scheduling, containerised tools), not a prediction-behaviour change. (`87b9bdc`, `af9fb67`, `9b43d59`)

## Benchmark (A. thaliana, brassicales_odb12, OMArk-enabled)

- gffcompare vs Araport11: Locus Sn/Pr **70.2 / 63.5**, Exon 80.5 / 84.1, Base 94.5 / 83.9 — matches prior runs within noise.
- OMArk on Brassicaceae HOGs (n = 17,996): Single **91.10 %**, Duplicated 7.85 % (Expected 5.51 %, Unexpected 2.33 %).

## Full commit list since v0.3.0-beta

```
026a390 fix: AUGUSTUS split-file rename missed genome_masked.fa stem
751c6ad add requirements.txt for host-side Python deps
166871e serialize OMArk via flock to avoid ete3 taxonomy SQLite contention
9b43d59 README benchmark: fix inconsistent bolding in comparison table
af9fb67 README benchmark: add F1 column
87b9bdc README benchmark: add galba.pl vs GALBA2 comparison
f3910c2 benchmark: declare A. thaliana genome as pre-softmasked
379c11c trim README redundancy; drop misleading "EP mode" naming
9de7d33 fix OMArk wiring: read omamer_db from [OMARK] section + env var
ba322dc removing redundant broken busco plot
9cdecc8 making busco work for offline nodes
```
