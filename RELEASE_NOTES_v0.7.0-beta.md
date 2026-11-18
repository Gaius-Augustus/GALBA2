# GALBA2 v0.7.0-beta — Release Notes

Release date: 2026-11-18
Previous release: [v0.6.0-beta](https://github.com/Gaius-Augustus/GALBA2/releases/tag/v0.6.0-beta)

Gene predictions are unchanged by this release: the toy scenario gives the
same gene set as earlier releases (258 genes, 0 differing lines in
`galba.gtf` and `galba.aa`). Most changes are ported from BRAKER4 (up to
BRAKER4 commit 7a14921, v0.8.0-beta).

## New features

- **GFF3 for ENA submission.** Set `ena_locus_tag_prefix` (and optionally
  `ena_taxid` and `ena_seqid_map`) in `samples.csv` to get
  `results/ena/galba.ena.gff3` for ENA decoupled annotation submission with
  Webin-CLI, plus an ID map and a report of every feature that was left out
  or changed ([scripts/gff3_to_ena.py](scripts/gff3_to_ena.py), rule
  `export_ena_gff3`, host Python, no container). The export keeps only INSDC
  qualifiers, assigns locus tags in genome order, includes the stop codon in
  the CDS, marks partial models and leaves out tRNAs and rRNAs that ENA
  rejects. See the README section "ENA GFF3 export for decoupled annotation
  submission". The prefix is checked at startup.
- **FANTASIA-Lite on CPU.** `[fantasia] device = cpu` (or
  `GALBA2_FANTASIA_DEVICE=cpu`) runs ProtT5 on `cpus_per_task` cores without
  a GPU request and without `fantasia_gpu_check`. Much slower than on a GPU;
  point `partition` at a CPU partition and allow a longer `max_runtime`.
  `device = gpu` stays the default.
- **Gzipped protein files.** `protein_fasta` may be `.gz`; such files are
  decompressed by `merge_proteins` before miniprot and the DIAMOND filter.
- **rRNA prediction with pybarrnap 0.5.1** instead of barrnap 0.9, which
  produced overlapping eukaryotic 5.8S/28S predictions (tseemann/barrnap#47).
  Container key `pybarrnap_image`.
- `galba.longest.gtf` and `galba.longest.aa` (longest coding isoform per
  gene) are copied to `results/` and gzipped.

## Bug fixes

- **DIAMOND filter emptied `galba.gtf`** when DIAMOND failed, for example on
  OrthoDB `*_train.fa` files that start with a blank line ("First line must
  begin with '>'"). `filter_gtf_by_diamond_against_ref.py` 1.2.0 runs DIAMOND
  on a cleaned copy of the proteins, stops when DIAMOND fails or no
  transcript has a hit, and keeps its temporary files next to the output.
  `sanity_check_augustus` fails when `galba.gtf` has no genes.
  `merge_proteins` drops blank lines and adds missing final newlines.
- **OOM kills of single-threaded rules** on nodes with many cores: they asked
  for `mem_of_node // cpus_per_task` (805 MB on a 256-CPU node). They now
  request at least 8 GB (never more than `mem_of_node`), more on each
  restart attempt (BRAKER4 eeae3cf).
- Default `galba_tools_image` is `galba2-tools:v0.5.0`; `galba-notebook:latest`
  has no `/usr/local/bin/miniprot`. `miniprot_index` stops with a clear
  message when the binary is missing.
- **Genome FASTA:** lines starting with `;` before the first sequence
  (FASTA comments, e.g. FAIR-bioHeaders metadata) are dropped; samtools
  faidx failed on them (Gaius-Augustus/BRAKER#982). A sequence name that
  occurs twice after cutting headers at the first space stops the run at
  `prepare_genome` instead of hours later (BRAKER#977).
- **Masking** (`run_masking`, `run_red_masking`) works in `[paths] tmp_dir`
  when it is set instead of always in `/tmp`, which can be too small for
  RepeatModeler2 on large genomes; the free space is logged.
- **`TMPDIR` not writable inside the job** (e.g. an SGE job directory that
  is not bound into the container): it is pointed at the job's scratch or
  fallback directory, so GNU sort does not fail on large inputs
  (BRAKER#861).
- **`run_augustus_hints`** stops with the AUGUSTUS error output of up to 5
  jobs when AUGUSTUS predicted no genes, instead of failing later without a
  clear message (BRAKER#887, #936).
- **`fix_in_frame_stop_codons`** passes the thread count to the script, which
  runs AUGUSTUS sequentially with one core or when Python multiprocessing is
  not available (no `/dev/shm`, WSL2; BRAKER#592, #873).
- **`normalize_gtf.py`** reads start and stop codons split by an intron and
  writes one codon line per part (BRAKER#751).
- **`create_new_species`** writes `logs/<sample>/create_new_species/create_new_species.log`
  and stops with the path when `new_species.pl` fails or writes no
  parameters file.

## Documentation

- README: x86_64 Linux only (no Apple Silicon); single-allocation mode for
  clusters where Snakemake cannot submit jobs from the login or compute
  nodes; gzipped proteins; `TMPDIR` handling; masking scratch; ENA export;
  FANTASIA on GPU or CPU (also in the overview figure).

## Known limitations

- Tested locally on the toy data (end to end, including gzipped proteins,
  a genome with `;` comment lines and the ENA export) and with unit tests
  (91 tests). Not re-run on a cluster for this release.
- FANTASIA with `device = cpu` was checked with Snakemake dry runs only.
- The ENA GFF3 has not been validated with Webin-CLI against a real
  submission; run `webin-cli -validate` first.
