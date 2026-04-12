# Migrating from GALBA (`galba.pl`) to GALBA2

This is a step-by-step guide for users who already know how to run the original GALBA Perl pipeline (`galba.pl`) and want to switch to GALBA2 (the Snakemake rewrite). It assumes you can already run a successful `galba.pl` annotation and that you have a working Singularity (or Apptainer) installation.

If you have never used GALBA before, read the main [README.md](README.md) first — this guide will not re-explain what GALBA does conceptually.

## What changes (and what does not)

| Concept | GALBA (`galba.pl`) | GALBA2 |
|---|---|---|
| Main entry point | `galba.pl --genome=... --prot_seq=...` | `snakemake --snakefile path/to/Snakefile ...` |
| Input specification | Command-line flags | Two text files: `samples.csv` and `config.ini` |
| Multiple genomes | One run per genome | Many genomes in one `samples.csv`, all run from one `snakemake` invocation |
| Output location | `--workingdir` | `output/{sample_name}/results/` |
| Resume after failure | Not possible — restart from scratch | Re-run the same `snakemake` command — picks up where it left off |
| HPC submission | You wrap `galba.pl` in your own SLURM script | `snakemake --executor slurm` submits each rule as a separate cluster job |
| Container | One container, you call it manually | All containers managed by Snakemake's `--use-singularity` |
| Tool paths | You configure `--AUGUSTUS_BIN_PATH`, `--MINIPROT_PATH`, etc. | Not needed — everything lives in the containers |
| QC outputs | None integrated | BUSCO, compleasm, OMArk, gffcompare integrated; HTML report generated automatically |

What does **not** change:

- The prediction logic. Miniprot aligns proteins, miniprothint extracts training genes and hints, AUGUSTUS is trained on the protein-derived genes and predicts with protein hints. The same algorithms with the same parameters.
- Output formats. The final files are still GTF, GFF3, protein FASTA, and coding sequences — just collected into the `results/` subdirectory and gzipped.

## Step 0 — prerequisites

Make sure you have:

- **Snakemake ≥ 7** (`snakemake --version`). We recommend Snakemake 8.x for SLURM support.
- **Singularity (or Apptainer) 3.x**. The Snakemake `--use-singularity` flag will pull all container images automatically on first run.
- **`snakemake-executor-plugin-slurm`** if you want to run on a SLURM cluster.
- **About 3.6 GB of disk for the Singularity image cache.** The main GALBA container is ~2 GB, plus BUSCO (~800 MB), AGAT (~370 MB), and smaller helpers. Plus enough disk for your genome FASTAs, intermediates, and outputs.

## Step 1 — clone the repository

```bash
git clone https://github.com/Gaius-Augustus/GALBA2.git
cd GALBA2
```

You will run snakemake from a working directory of your choice and point it at the `Snakefile` here. You do not need to install GALBA2 — there is nothing to install. The Snakefile is the entry point.

## Step 2 — translate your `galba.pl` command

Take the command you currently use, and write down each flag. Then map each flag to either a `samples.csv` column or a `config.ini` parameter using the table below.

### Common `galba.pl` flag → GALBA2 destination

| `galba.pl` flag | GALBA2 destination |
|---|---|
| `--genome=genome.fa` | `samples.csv` column `genome` |
| `--prot_seq=proteins.fa` | `samples.csv` column `protein_fasta` (colon-separated for multiple files) |
| `--species=sname` | `samples.csv` column `sample_name` |
| `--workingdir=/path/to/wd/` | The directory you `cd` into before running `snakemake` (output goes to `output/{sample_name}/`) |
| `--threads=N` | `snakemake --cores N` (local) or `--cores N --jobs N` (SLURM) |
| `--gff3` | Always on. Both `.gtf` and `.gff3` are produced. |
| `--skipOptimize` | `config.ini` `[PARAMS]` `skip_optimize_augustus = 1` |
| `--nocleanup` | `config.ini` `[PARAMS]` `no_cleanup = 1` |
| `--translation_table=N` | `config.ini` `[PARAMS]` `translation_table = N` |
| `--disable_diamond_filter` | `config.ini` `[PARAMS]` `disable_diamond_filter = 1` |
| `--splice_sites=gcag,atac` | `config.ini` `[PARAMS]` `allow_hinted_splicesites = gcag,atac` |
| `--eval=reference.gtf` | `samples.csv` column `reference_gtf` |
| `--softmasking` | Default. No flag needed. |
| `--AUGUSTUS_CONFIG_PATH=...` | `config.ini` `[paths]` `augustus_config_path = ...` (default `augustus_config` works for most users) |
| `--AUGUSTUS_BIN_PATH`, `--MINIPROT_PATH`, `--DIAMOND_PATH`, etc. | Not needed. Everything lives in the containers. |

### Flags that do not have a GALBA2 equivalent

Some `galba.pl` options were not ported to GALBA2. This is intentional — they were either rarely used, superseded by better defaults, or no longer relevant in a containerized Snakemake workflow:

- `--prg=gth` — GenomeThreader is no longer supported. GALBA2 uses miniprot exclusively.
- `--AUGUSTUS_ab_initio` — Ab initio predictions are not produced separately.
- `--crf` / `--keepCrf` — CRF training is not supported; GALBA2 uses HMM training only.
- `--makehub` / `--email` — UCSC track hub generation was removed.
- `--hints=file.gff` — Pre-generated hints input is not supported; GALBA2 generates hints from protein alignments.
- `--traingenes=file.gtf` — Pre-generated training genes are not supported.
- `--augustus_args="..."` — Custom AUGUSTUS arguments are not exposed.
- `--rounds=N` — Optimization rounds are handled internally.
- `--checkSoftware` — Use `snakemake --dryrun` instead.
- `--verbosity=N` — Use Snakemake's own `--quiet` or `-v` flags.

If you relied on any of these, the main [README.md](README.md) explains the GALBA2 alternatives where applicable.

## Step 3 — create your `samples.csv`

Pick a working directory (it does not have to be inside the cloned repo):

```bash
mkdir ~/my_galba2_run
cd ~/my_galba2_run
```

Write a `samples.csv`. The header is fixed — copy it exactly:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
```

Then add **one row per genome**. Leave a column empty if you do not have that input.

### Typical example

If you used to run:

```bash
galba.pl --genome=genome.fa --prot_seq=proteins.fa \
         --species=my_species --threads=8
```

Your `samples.csv` is:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,/path/to/proteins.fa,eukaryota_odb12,
```

Notes:

- **`busco_lineage` is required.** Pick a clade-specific lineage like `arthropoda_odb12`, `viridiplantae_odb12`, etc. If you want the most generic option, use `eukaryota_odb12`.
- **Paths can be absolute** (recommended) **or relative** to the working directory you run `snakemake` from.
- The `sample_name` becomes the AUGUSTUS species name and the output directory name. Pick something stable — re-running with the same name resumes; changing it forces a fresh run.

### With a pre-masked genome

If you used to provide a masked genome:

```bash
galba.pl --genome=genome.fa --prot_seq=proteins.fa --softmasking \
         --species=my_species --threads=8
```

Put the masked genome in the `genome_masked` column:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,/path/to/genome_masked.fa,/path/to/proteins.fa,eukaryota_odb12,
```

### With a reference annotation for evaluation

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,/path/to/proteins.fa,arthropoda_odb12,/path/to/reference.gtf
```

### Multi-sample example

You can annotate several genomes in the same invocation — something `galba.pl` cannot do:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
fly,/data/fly.fa,,/data/fly_proteins.fa,arthropoda_odb12,
plant,/data/plant.fa,/data/plant_masked.fa,/data/plant_proteins.fa,viridiplantae_odb12,/data/plant_ref.gtf
```

## Step 4 — create your `config.ini`

In the same working directory, create `config.ini`:

```ini
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[PARAMS]
skip_optimize_augustus = 0
disable_diamond_filter = 0
skip_busco = 0
run_omark = 0
translation_table = 1
no_cleanup = 0

[SLURM_ARGS]
cpus_per_task = 48
mem_of_node = 120000
max_runtime = 4320
```

The defaults above are appropriate for most runs on a 48-core SLURM node with 120 GB RAM and a 72-hour wall time. Adjust to your cluster.

A few flags you might want to flip:

| Flag | When to set it |
|---|---|
| `skip_optimize_augustus = 1` | Skip the slow AUGUSTUS parameter optimization step. Saves time but slightly reduces accuracy. Recommended only for testing. |
| `skip_busco = 1` | Skip BUSCO completeness assessment. Saves ~20 min. |
| `run_omark = 1` | Run OMArk on the predicted proteins. Requires the LUCA.h5 database (~15 GB). |
| `disable_diamond_filter = 1` | Skip DIAMOND filtering of predictions against input proteins. |
| `no_cleanup = 1` | Keep all intermediate files. Useful for debugging but disk-hungry. |

Every value in `[PARAMS]` and `[SLURM_ARGS]` can also be overridden via an environment variable named `GALBA2_<KEY_UPPER>`, e.g.:

```bash
GALBA2_SKIP_OPTIMIZE_AUGUSTUS=1 GALBA2_MAX_RUNTIME=240 snakemake ...
```

For the full list see the README's [Description of selected configuration options](README.md#description-of-selected-configuration-options).

## Step 5 — run GALBA2

Two cases — local workstation or HPC cluster.

### Local (workstation, no SLURM)

From your working directory (the one with `samples.csv` and `config.ini`):

```bash
snakemake \
    --snakefile /path/to/GALBA2/Snakefile \
    --cores 8 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home"
```

- `--cores 8` → use 8 CPU threads.
- `--use-singularity` → run all rules inside their declared containers.
- `--singularity-prefix .singularity_cache` → cache pulled images here. Without this flag, every working directory gets its own copy of the container images.
- `--singularity-args "-B /home"` → make `/home` visible inside the container. **If your data lives outside `/home`** (e.g. `/scratch`, `/data`, `/gpfs`), add those bind paths: `--singularity-args "-B /home -B /scratch -B /data"`. The container can only see directories that are explicitly bound.

The first run will pull all containers GALBA2 needs (~3.6 GB total). Give it a few minutes depending on your network. Subsequent runs reuse the cached images.

### HPC cluster (SLURM)

From your working directory:

```bash
snakemake \
    --snakefile /path/to/GALBA2/Snakefile \
    --executor slurm \
    --default-resources slurm_partition=batch mem_mb=120000 \
    --cores 48 \
    --jobs 48 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home -B /scratch"
```

- `--executor slurm` → submit each rule as a separate `sbatch` job. Set `slurm_partition=` to your cluster's partition name.
- `--cores 48` → cores per individual SLURM job. Match your `cpus_per_task` in `config.ini`.
- `--jobs 48` → maximum number of SLURM jobs in flight at once.
- `mem_mb=120000` → memory per SLURM job (MB). Match `mem_of_node` in `config.ini`.

**Run snakemake itself in `screen` or `tmux`**, or submit it as a long-running SLURM job. The Snakemake driver process must stay alive for the entire pipeline duration.

### Resume after a failure

If anything fails (out of memory, time limit, network glitch), just **re-run the same command**. Snakemake reads the existing partial outputs, figures out what is missing, and only re-runs what is needed. There is no `--useexisting` flag — resume is the default behavior.

## Step 6 — find your outputs

Your final files are in:

```
output/{sample_name}/results/
├── galba.gtf.gz          # final gene set, GTF
├── galba.gff3.gz         # same set, GFF3 (AGAT-formatted)
├── galba.aa.gz           # protein sequences
├── galba.codingseq.gz    # coding sequences
├── hintsfile.gff.gz      # protein-derived hints
├── galba_report.html     # self-contained HTML report
├── galba_citations.bib   # BibTeX for all tools used
├── software_versions.tsv # versions of all software used
└── quality_control/
    ├── busco_summary.txt
    ├── compleasm_summary.txt
    ├── completeness.png
    ├── gene_set_statistics.txt
    └── (more plots and tables)
```

The main differences from `galba.pl` output:

1. **Outputs are gzipped** by default in `results/`.
2. **Outputs are in a `results/` subdirectory.** Everything else in `output/{sample_name}/` is intermediate state.
3. **The HTML report is new.** Open `galba_report.html` in a browser — it embeds all QC plots, gene-set statistics, software versions, and citations into a single self-contained file.
4. **The BibTeX file is new.** All software cited in your run, ready to drop into a bibliography manager.

## Step 7 — sanity-check the run

Open `galba_report.html` in a browser. It will tell you:

- How many genes were predicted at each pipeline stage (training set → AUGUSTUS → final).
- BUSCO and compleasm completeness scores against your `busco_lineage`.
- Gene-set statistics (CDS lengths, exons per gene, mono-exonic vs multi-exonic ratio, etc.).
- The exact software versions and citations.

If a step looks wrong, look at the per-rule logs in `logs/{sample_name}/<rule_name>/` to find what happened. On SLURM, the captured stdout/stderr of each job lives in `.snakemake/slurm_logs/rule_<rule_name>/<sample_name>/<jobid>.log`.

## Common gotchas when migrating

### "Singularity cannot find my data"

By default `--singularity-args "-B /home"` only binds `/home`. If your genome FASTA is in `/scratch/myproject/genome.fa`, the container cannot see it. Add the bind path:

```
--singularity-args "-B /home -B /scratch"
```

This is the most common first-run failure.

### "I don't have a `busco_lineage` for my organism"

Use the most specific BUSCO odb12 lineage that contains your organism. If you really do not know which to pick, `eukaryota_odb12` works for any eukaryote (it is just less informative for QC than a clade-specific one).

### "I want to override one config value without editing config.ini"

Every value in `[PARAMS]` and `[SLURM_ARGS]` can be overridden via an environment variable named `GALBA2_<KEY_UPPER>`. For example:

```bash
GALBA2_SKIP_OPTIMIZE_AUGUSTUS=1 GALBA2_MAX_RUNTIME=240 snakemake ...
```

### "Snakemake says my SLURM job has no `mem_mb` resource"

Add `mem_mb=...` to `--default-resources`:

```
--default-resources slurm_partition=batch mem_mb=120000
```

`galba.pl` did not need this because you set memory yourself in the wrapping `sbatch` script. In GALBA2, Snakemake submits each rule individually and needs to know how much memory to request.

### "I want to keep all intermediate files for debugging"

Set `no_cleanup = 1` in `[PARAMS]`. By default the `collect_results` rule moves the important outputs into `results/` and deletes the intermediate state. With `no_cleanup = 1`, everything in `output/{sample_name}/` is preserved.

### "I had a `--workingdir` separate from my source code"

In GALBA2, your "working directory" is whatever directory you `cd` into before running snakemake. It contains `samples.csv` and `config.ini`, and Snakemake creates `output/`, `logs/`, `benchmarks/`, `.snakemake/` underneath it. Pick any directory you like — it does not have to be inside the cloned GALBA2 repository. Reference the Snakefile via its absolute path (`--snakefile /path/to/GALBA2/Snakefile`).

### "My runs hung / SLURM showed COMPLETING jobs forever"

If you killed a snakemake driver mid-run, leftover SLURM jobs may still be in the queue and stale lock files may be in `.snakemake/locks/`. Clean both before re-running:

```bash
# scancel any stale jobs (use squeue -u $USER to identify)
scancel <jobid> ...
# remove stale locks
rm -f .snakemake/locks/*
# now re-run snakemake
snakemake ...
```

## Where to read more

- [README.md](README.md) — full reference: every config option, every output file.
- [test_scenarios/](test_scenarios/) — runnable HPC test scenarios.
- [test_scenarios_local/](test_scenarios_local/) — same scenarios but for a local workstation (no SLURM).

If you get stuck, the per-rule logs in `logs/{sample_name}/` and the SLURM job logs in `.snakemake/slurm_logs/` are usually enough to identify what went wrong. Open an issue at <https://github.com/Gaius-Augustus/GALBA2/issues> if you cannot.
