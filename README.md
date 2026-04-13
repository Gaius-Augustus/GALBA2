<p align="center"><img src="img/logo.png" width="280" height="280" alt="GALBA2 logo"></p>

<p align="center">
  <a href="https://github.com/Gaius-Augustus/GALBA2/releases"><img src="https://img.shields.io/github/v/release/Gaius-Augustus/GALBA2?include_prereleases&sort=semver" alt="Release"></a>
  <a href="LICENSE"><img src="https://img.shields.io/github/license/Gaius-Augustus/GALBA2" alt="License"></a>
  <img src="https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen" alt="Snakemake &ge; 7.0">
  <img src="https://img.shields.io/badge/container-singularity-orange" alt="Singularity">
</p>

# GALBA2 — protein homology based genome annotation in Snakemake

Author: Katharina J. Hoff

Contact for Repository
========================

Katharina J. Hoff, University of Greifswald, Germany, katharina.hoff@uni-greifswald.de, +49 3834 420 4624

> **Migrating from GALBA (`galba.pl`)?** The fastest path to your first GALBA2 run is the dedicated step-by-step tutorial: **[MIGRATING_FROM_GALBA.md](MIGRATING_FROM_GALBA.md)**. It maps every `galba.pl` flag you already use to the equivalent `samples.csv` column or `config.ini` parameter, and walks through the translation with concrete examples.

Contents
========

-   [What is different in GALBA2?](#what-is-different-in-galba2)
-   [What is GALBA?](#what-is-galba)
-   [Keys to successful gene prediction](#keys-to-successful-gene-prediction)
-   [Protein database preparation](#protein-database-preparation)
-   [Installation](#installation)
    -   [Snakemake](#snakemake)
    -   [Singularity](#singularity)
-   [Running GALBA2](#running-galba2)
    -   [Preparing input files](#preparing-input-files)
        -   [samples.csv](#samplescsv)
        -   [config.ini](#configini)
    -   [Running locally](#running-locally)
    -   [Running on an HPC cluster with SLURM](#running-on-an-hpc-cluster-with-slurm)
    -   [Description of selected configuration options](#description-of-selected-configuration-options)
-   [Output of GALBA2](#output-of-galba2)
-   [Example data](#example-data)
-   [Bug reporting](#bug-reporting)
-   [Citing GALBA2 and software called by GALBA2](#citing-galba2-and-software-called-by-galba2)
-   [License](#license)

What is different in GALBA2?
=============================

GALBA2 is a rewrite of the [GALBA pipeline](https://github.com/Gaius-Augustus/GALBA) in Snakemake, based on the [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) framework. The gene prediction logic is the same: miniprot aligns proteins to the genome, miniprothint extracts training genes and protein hints, AUGUSTUS is trained on the protein-derived genes and predicts with protein hints. What changed is how this logic is orchestrated.

**The old GALBA** was a monolithic Perl script (`galba.pl`, ~275 KB). It managed all tool calls, error handling, and file paths in a single script. It could not resume after failures and did not natively support cluster execution.

**GALBA2** replaces that script with a Snakemake workflow. All bioinformatics tools run inside a Singularity container. You do not install miniprot, AUGUSTUS, DIAMOND, or any of their dependencies on your system. Snakemake handles job scheduling, parallelization, and automatic resume after failures.

Key differences:

-   **CSV-based multi-sample input.** You can annotate multiple genomes in a single run. Each row in `samples.csv` defines one genome with its protein evidence.

-   **Automatic resume.** If a run fails (out of memory, time limit, network error), re-run the same command. Snakemake picks up where it left off.

-   **Modular rules.** Each step is a separate Snakemake rule in its own file. This makes it straightforward to understand, debug, and extend the pipeline.

-   **HPC-ready with SLURM.** Snakemake's SLURM executor submits each rule as a separate cluster job.

-   **Integrated postprocessing and quality control.** GFF3 conversion (AGAT), DIAMOND filtering against input proteins, BUSCO/compleasm completeness assessment, and optional evaluation against a reference annotation (gffcompare) are included.

Benchmark accuracy vs native galba.pl
======================================

To verify that GALBA2 reproduces the gene-prediction accuracy of the original `galba.pl` Perl pipeline, we run the same *Arabidopsis thaliana* genome (TAIR10 assembly, ~121 Mb) through both pipelines with matched configurations and score the resulting gene sets against the Phytozome Araport11 reference annotation using `gffcompare` v0.12.6 at the CDS level (`--strict-match -e 3 -T`).

**Inputs (identical for all pipelines):** TAIR10 genome FASTA (pre-softmasked), Viridiplantae proteins from OrthoDB v12 (~467 MB). All pipelines run with `optimize_augustus.pl` enabled. BRAKER4 EP (which uses ProtHint + GeneMark-EP+ instead of miniprot) is included as an additional reference point.

<table>
  <thead>
    <tr>
      <th colspan="2" align="center">Locus Sn</th>
      <th colspan="2" align="center">Locus Pr</th>
      <th colspan="2" align="center">Exon Sn</th>
      <th colspan="2" align="center">Exon Pr</th>
      <th colspan="2" align="center">Base Sn</th>
      <th colspan="2" align="center">Base Pr</th>
    </tr>
    <tr>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
      <th>GALBA2</th><th>BRAKER4&nbsp;EP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="right"><b>70.4</b></td><td align="right">77.0</td>
      <td align="right"><b>63.4</b></td><td align="right">61.8</td>
      <td align="right">80.8</td><td align="right"><b>83.2</b></td>
      <td align="right"><b>84.3</b></td><td align="right">80.4</td>
      <td align="right">94.7</td><td align="right"><b>96.1</b></td>
      <td align="right"><b>84.1</b></td><td align="right">80.4</td>
    </tr>
  </tbody>
</table>

Sn = Sensitivity (% of reference features recovered). Pr = Precision (% of predicted features matching reference). Higher is better for both. Bold marks the better value per metric.

**Reading the table.** GALBA2 predicts 30,747 genes with higher precision than BRAKER4 EP (34,834 genes) across all levels: locus precision 63.4 vs 61.8, exon precision 84.3 vs 80.4, base precision 84.1 vs 80.4. BRAKER4 EP achieves higher sensitivity because ProtHint + GeneMark-EP+ cast a wider net, but at the cost of more false positives. When closely related protein evidence is available, GALBA2 provides a favourable balance of sensitivity and precision.

**Where GALBA shines.** The *A. thaliana* genome (~121 Mb) is small. On larger genomes — which are common for many animal, plant, and fungal species — GALBA's miniprot-based approach scales much better than BRAKER2's ProtHint + GeneMark-EP+ pipeline. On genomes of 1 Gb and above, GALBA consistently outperforms BRAKER2 EP mode in both accuracy and runtime, because miniprot handles large, repeat-rich genomes more robustly than the Spaln/DIAMOND-based ProtHint alignment. If you are annotating a large genome with protein evidence only, GALBA2 is the recommended tool.

What is GALBA?
===============

GALBA is a pipeline for genome annotation using protein homology. It is designed for cases where no RNA-Seq data is available, but a database of protein sequences from related species exists. GALBA uses miniprot to align proteins to the genome, generates training genes from high-quality alignments, trains AUGUSTUS on these genes, and predicts genes using protein-derived hints.

GALBA is particularly useful for:

-   **Novel species** without RNA-Seq data, where a protein database from related species is available.
-   **Large-scale annotation projects** where dozens of genomes need to be annotated with the same protein database.
-   **Species with high-quality protein databases** (e.g. from OrthoDB) but no transcriptome evidence.

<p align="center">
  <img src="img/pipeline_overview.svg" alt="GALBA2 pipeline overview" width="85%">
  <br>
  <em>Figure&nbsp;1. Overview of the GALBA2 pipeline. Solid borders mark required steps; dashed borders mark optional steps. Miniprot aligns proteins to the genome, miniprothint extracts training genes and hints, AUGUSTUS is trained and predicts genes with protein evidence.</em>
</p>

Keys to successful gene prediction
====================================

-   Use a high quality genome assembly. If you have a huge number of very short scaffolds, those will increase runtime but not prediction accuracy.

-   Use simple scaffold names in the genome file (e.g. `>contig1`). Complex headers can cause parsing issues.

-   The genome should be soft-masked for repeats. Soft-masking (lowercase letters for repeats) leads to better results than hard-masking (replacing with `N`). If you provide an unmasked genome, GALBA2 will still run but prediction quality may suffer in repetitive regions.

-   Always check gene prediction results before further usage. You can use a genome browser for visual inspection of gene models in context with protein alignments.

Protein database preparation
=============================

GALBA2 requires protein sequences from **closely related species** — this is different from BRAKER, which uses large OrthoDB partitions with distant homologs. GALBA relies on miniprot to produce high-quality spliced alignments, which works best with closely related proteins that share high sequence identity with the target genome.

**Recommended protein sources:**

-   Annotated proteomes from **4 or more closely related species** (e.g. from the same family or order). The more species, the better the coverage, but all should be reasonably close relatives.
-   Each protein must have a **unique identifier** in the FASTA header. If you combine proteomes from multiple species, ensure there are no duplicate IDs.
-   **Do NOT use large OrthoDB partitions** (e.g. all of Viridiplantae or Metazoa). These contain proteins from very distant species that produce poor miniprot alignments and lead to unreliable training genes. For such cases, use [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) instead, which uses ProtHint/GeneMark and handles distant homologs.

**Example:** To annotate a new grass species, combine the proteomes of *Oryza sativa*, *Zea mays*, *Sorghum bicolor*, and *Brachypodium distachyon* into a single FASTA file:

```bash
cat rice_proteins.fa maize_proteins.fa sorghum_proteins.fa brachy_proteins.fa > close_relatives.fa
```

Then specify it in `samples.csv`:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
my_grass,genome.fa,,close_relatives.fa,poales_odb12,
```

You can also provide multiple files colon-separated:

```csv
protein_fasta
rice_proteins.fa:maize_proteins.fa:sorghum_proteins.fa
```

Installation
============

**Platform requirement:** GALBA2 requires **Linux**. It relies on GNU coreutils (`readlink -f`), GNU grep (`grep -P`), and Singularity/Apptainer. If you are on macOS or Windows, run GALBA2 inside a Linux virtual machine or WSL2.

GALBA2 requires two things on your system: Snakemake and Singularity. All bioinformatics tools (miniprot, miniprothint, AUGUSTUS, DIAMOND, etc.) run inside the container. You do not need to install them.

Snakemake
---------

We recommend installing Snakemake with `pip` into a virtual environment:

```
python3 -m venv snakemake_env
source snakemake_env/bin/activate
pip install snakemake
```

If you intend to run GALBA2 on an HPC cluster with SLURM, you also need the Snakemake SLURM executor plugin:

```
pip install snakemake-executor-plugin-slurm
```

Singularity
-----------

Singularity (or Apptainer, its open-source successor) must be installed on your system. On most HPC clusters, it is already available as a module:

```
module load singularity
```

GALBA2 will automatically pull the container image and convert it to a Singularity image on the first run.

**Container image:**

| Container | Image | Size | Used for |
| --- | --- | --- | --- |
| GALBA2 tools | `katharinahoff/galba2-tools:latest` | 112 MB | miniprot, miniprot_boundary_scorer, miniprothint, compleasm, matplotlib, BioPython — protein alignment, hint generation, completeness assessment, report plots |
| AUGUSTUS | `quay.io/biocontainers/augustus:3.5.0--pl5321h9716f88_9` | 249 MB | AUGUSTUS, etraining, DIAMOND, all Perl training/prediction scripts — training, prediction, postprocessing |
| AGAT | `quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0` | 370 MB | GTF-to-GFF3 conversion |
| TE Tools | `dfam/tetools:latest` | 727 MB | RepeatModeler2 + RepeatMasker (only when genome masking is needed) |
| BUSCO | `ezlabgva/busco:v6.0.0_cv1` | 801 MB | BUSCO completeness assessment (optional) |
| OMArk | `quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0` | 455 MB | OMArk proteome quality assessment (optional, only when `run_omark = 1`) |
| gffcompare | `quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1` | 11 MB | Evaluation against a reference annotation (optional) |
| FANTASIA-Lite | `katharinahoff/fantasia_for_brain:lite.v0.0.2` | ~6 GB | Functional GO annotation via ProtT5 embeddings (optional, **GPU-only**) |

A standard GALBA2 run (no masking, no BUSCO) pulls only the GALBA2 tools + AUGUSTUS + AGAT containers: **~730 MB** total.

**Singularity bind paths:** Singularity containers can only see directories that are explicitly bound. GALBA2 passes `--singularity-args "-B /home"` by default. If your data resides outside `/home` (e.g. in `/scratch`, `/data`, or `/gpfs`), you must add those paths:

```
--singularity-args "-B /home -B /scratch -B /data"
```

Running GALBA2
===============

Preparing input files
---------------------

GALBA2 requires two configuration files in your working directory: `samples.csv` and `config.ini`.

### samples.csv

This CSV file defines your input samples. Each row is one genome to annotate.

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
```

**Columns:**

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | yes | Unique identifier. Output goes to `output/{sample_name}/`. |
| `genome` | yes | Path to genome FASTA file. |
| `genome_masked` | no | Path to soft-masked genome. If empty, the unmasked genome is used directly. |
| `protein_fasta` | **yes** | Protein sequences in FASTA format. Multiple files can be colon-separated. We recommend OrthoDB. |
| `busco_lineage` | **yes** | BUSCO lineage for QC (e.g. `eukaryota_odb12`, `arthropoda_odb12`). |
| `reference_gtf` | no | Reference annotation for gffcompare evaluation. Optional. |

**Example: annotating two genomes in one run**

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
fly,fly_genome.fa,,orthodb_arthropoda.fa,arthropoda_odb12,
plant,plant_genome.fa,plant_masked.fa,orthodb_viridiplantae.fa,viridiplantae_odb12,reference.gtf
```

You can also combine multiple protein sources by colon-separating them:

```csv
sample_name,genome,genome_masked,protein_fasta,busco_lineage,reference_gtf
my_species,genome.fa,,orthodb.fa:close_relative.fa,eukaryota_odb12,
```

### config.ini

This file contains pipeline parameters. Place it in the same directory as your `samples.csv`.

```ini
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config
; busco_download_path = shared_data/busco_downloads      # optional, auto-discovered
; compleasm_download_path = shared_data/compleasm_downloads  # optional, auto-discovered

[containers]
galba_tools_image = docker://katharinahoff/galba2-tools:latest    # miniprot, miniprothint, compleasm
augustus_image = docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9716f88_9  # AUGUSTUS, DIAMOND

[PARAMS]
skip_optimize_augustus = 0          # set to 1 to skip AUGUSTUS optimization (saves time)
disable_diamond_filter = 0         # set to 1 to skip DIAMOND filtering of predictions
skip_busco = 0                     # set to 1 to skip the BUSCO pipeline
run_omark = 0                      # set to 1 to run OMArk (requires LUCA.h5 database)
translation_table = 1              # genetic code table (1=standard)
allow_hinted_splicesites = gcag,atac # non-canonical splice sites for AUGUSTUS
augustus_chunksize = 3000000        # genome chunk size (bp) for parallel AUGUSTUS prediction
augustus_overlap = 500000           # overlap (bp) between adjacent AUGUSTUS chunks
no_cleanup = 0                     # set to 1 to keep all intermediate files (debugging only)
use_dev_shm = 0                    # set to 1 to use /dev/shm for temp files (faster I/O)

[SLURM_ARGS]
cpus_per_task = 48
mem_of_node = 120000                # memory in MB
max_runtime = 4320                  # runtime in minutes (72 hours)
```

Every value in `[PARAMS]` and `[SLURM_ARGS]` can also be overridden via an environment variable named `GALBA2_<KEY_UPPER>`, e.g. `GALBA2_MAX_RUNTIME=240` or `GALBA2_SKIP_OPTIMIZE_AUGUSTUS=1`. Environment variables win over the file. The path of the file itself can be overridden with `GALBA2_CONFIG=/path/to/another.ini`.

Running locally
---------------

For running GALBA2 on a local workstation or a single compute node:

```
snakemake --cores 8 --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home"
```

Adjust `--cores` to the number of CPU cores available.

Running on an HPC cluster with SLURM
--------------------------------------

For running on a SLURM-managed cluster:

```
snakemake \
    --executor slurm \
    --default-resources slurm_partition=batch mem_mb=120000 \
    --cores 48 \
    --jobs 48 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home -B /scratch"
```

We recommend running the Snakemake master process in a `screen` or `tmux` session, or submitting it as a long-running SLURM job.

Description of selected configuration options
----------------------------------------------

| Option | Default | Description |
|--------|---------|-------------|
| `skip_optimize_augustus` | 0 | Skip `optimize_augustus.pl`. Saves hours on large genomes. **Do not enable for production runs** — optimization improves accuracy. |
| `disable_diamond_filter` | 0 | Skip DIAMOND filtering of AUGUSTUS predictions against input proteins. The filter removes gene predictions that have no sequence similarity match in the protein database. Disabling produces more genes but also more false positives. |
| `skip_busco` | 0 | Skip the BUSCO completeness assessment (slow). compleasm still runs. |
| `run_omark` | 0 | Run OMArk proteome quality assessment. Requires the LUCA.h5 OMAmer database (~8.8 GB). Download with `bash scripts/download_data.sh --omark`. |
| `run_fantasia` | 0 | Run FANTASIA-Lite functional GO annotation. **GPU-only** (requires NVIDIA GPU + pre-staged container and ProtT5 cache). See `[fantasia]` section in `config.ini`. |
| `translation_table` | 1 | NCBI genetic code table. Table 1 is the standard code. |
| `no_cleanup` | 0 | Keep all intermediate files. Useful for debugging. |
| `use_dev_shm` | 0 | Use `/dev/shm` for temporary files during AUGUSTUS prediction (faster I/O on nodes with large RAM). |

### Downloading shared databases

GALBA2 ships a helper script to download shared databases:

```bash
# Download compleasm/BUSCO lineage cache directories
bash scripts/download_data.sh

# Also download the OMAmer database for OMArk
bash scripts/download_data.sh --omark

# Use a custom location (shared with BRAKER4 or other pipelines)
bash scripts/download_data.sh --shared-data /data/shared_databases
```

The script auto-discovers existing shared_data directories in sibling BRAKER4 or BRAKER-as-snakemake repositories to avoid redundant downloads.

Output of GALBA2
=================

GALBA2 collects all important output files in `output/{sample}/results/`:

| File | Description |
|------|-------------|
| `galba.gtf.gz` | Gene predictions in GTF format |
| `galba.gff3.gz` | Gene predictions in GFF3 format |
| `galba.aa.gz` | Predicted protein sequences |
| `galba.codingseq.gz` | Predicted coding sequences |
| `hintsfile.gff.gz` | Protein-derived hints used for prediction |
| `galba_report.html` | Pipeline report with statistics |
| `galba_citations.bib` | BibTeX citations for tools used |
| `software_versions.tsv` | Versions of all tools used |
| `quality_control/` | BUSCO/compleasm completeness, gene set statistics, training summary |

Example data
=============

GALBA2 ships with a small toy dataset in `test_data/` (a ~1 MB genome and ~119 KB protein set from the original GALBA test suite). To run the toy test locally:

```bash
cd test_scenarios_local/scenario_01_ep
bash run_test.sh
```

For HPC benchmarks, see `test_scenarios/scenario_benchmark_athaliana_ep/` which annotates the full *A. thaliana* genome with Viridiplantae proteins from OrthoDB.

Bug reporting
==============

Please report bugs and feature requests at https://github.com/Gaius-Augustus/GALBA2/issues

Before reporting a bug, please check:

-   Is your genome file in valid FASTA format with simple headers?
-   Is your protein database a multi-FASTA file with many representative proteins per family?
-   Do you have sufficient disk space? GALBA2 needs several GB for intermediate files.
-   Check the Snakemake log in `.snakemake/log/` for error details.

Citing GALBA2 and software called by GALBA2
=============================================

If you use GALBA2, please cite:

-   Hoff, K. J. (2026). GALBA2: a Snakemake pipeline for protein homology based genome annotation. *In preparation.*
-   Hoff, K. J., Bruna, T., Lomsadze, A., Borodovsky, M., & Stanke, M. (2023). A pipeline for automated gene prediction for novel species using protein homology. *NAR Genomics and Bioinformatics*, 5(3), lqad064. doi:10.1093/nargab/lqad064

Software called by GALBA2:

-   Li, H. (2023). Protein-to-genome alignment with miniprot. *Bioinformatics*, 39(1), btad014. doi:10.1093/bioinformatics/btad014
-   Stanke, M., Diekhans, M., Baertsch, R., & Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. *Bioinformatics*, 24(5), 637-644. doi:10.1093/bioinformatics/btn013
-   Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods*, 18, 366-368. doi:10.1038/s41592-021-01101-x
-   Huang, N., & Li, H. (2023). compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics*, 39(10), btad595. doi:10.1093/bioinformatics/btad595
-   Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. doi:10.12688/f1000research.29032.2

Usage of AI
===========

The development of GALBA2 was assisted by AI tools:

-   **Claude** (Anthropic) was used for code generation and pipeline architecture. Most Snakemake rules, helper scripts, and the overall workflow design were developed in collaboration with Claude.
-   **Gemini** (Google) was used to generate the GALBA2 logo.

All AI-generated code was reviewed and tested by the authors. The pipeline was validated end-to-end on test data and the *A. thaliana* benchmark genome.

License
=======

See [LICENSE](LICENSE) for details.
