# Graph Report - /home/katharina/git/GALBA2  (2026-04-22)

## Corpus Check
- 19 files · ~528,323 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 294 nodes · 328 edges · 76 communities detected
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 12 edges (avg confidence: 0.83)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 28|Community 28]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 32|Community 32]]
- [[_COMMUNITY_Community 33|Community 33]]
- [[_COMMUNITY_Community 34|Community 34]]
- [[_COMMUNITY_Community 35|Community 35]]
- [[_COMMUNITY_Community 36|Community 36]]
- [[_COMMUNITY_Community 37|Community 37]]
- [[_COMMUNITY_Community 38|Community 38]]
- [[_COMMUNITY_Community 39|Community 39]]
- [[_COMMUNITY_Community 40|Community 40]]
- [[_COMMUNITY_Community 41|Community 41]]
- [[_COMMUNITY_Community 42|Community 42]]
- [[_COMMUNITY_Community 43|Community 43]]
- [[_COMMUNITY_Community 44|Community 44]]
- [[_COMMUNITY_Community 45|Community 45]]
- [[_COMMUNITY_Community 46|Community 46]]
- [[_COMMUNITY_Community 47|Community 47]]
- [[_COMMUNITY_Community 48|Community 48]]
- [[_COMMUNITY_Community 49|Community 49]]
- [[_COMMUNITY_Community 50|Community 50]]
- [[_COMMUNITY_Community 51|Community 51]]
- [[_COMMUNITY_Community 52|Community 52]]
- [[_COMMUNITY_Community 53|Community 53]]
- [[_COMMUNITY_Community 54|Community 54]]
- [[_COMMUNITY_Community 55|Community 55]]
- [[_COMMUNITY_Community 56|Community 56]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 58|Community 58]]
- [[_COMMUNITY_Community 59|Community 59]]
- [[_COMMUNITY_Community 60|Community 60]]
- [[_COMMUNITY_Community 61|Community 61]]
- [[_COMMUNITY_Community 62|Community 62]]
- [[_COMMUNITY_Community 63|Community 63]]
- [[_COMMUNITY_Community 64|Community 64]]
- [[_COMMUNITY_Community 65|Community 65]]
- [[_COMMUNITY_Community 66|Community 66]]
- [[_COMMUNITY_Community 67|Community 67]]
- [[_COMMUNITY_Community 68|Community 68]]
- [[_COMMUNITY_Community 69|Community 69]]
- [[_COMMUNITY_Community 70|Community 70]]
- [[_COMMUNITY_Community 71|Community 71]]
- [[_COMMUNITY_Community 72|Community 72]]
- [[_COMMUNITY_Community 73|Community 73]]
- [[_COMMUNITY_Community 74|Community 74]]
- [[_COMMUNITY_Community 75|Community 75]]

## God Nodes (most connected - your core abstractions)
1. `GALBA2` - 43 edges
2. `GALBA2 pipeline overview diagram` - 15 edges
3. `main()` - 13 edges
4. `main()` - 8 edges
5. `read_file()` - 7 edges
6. `generate_html()` - 7 edges
7. `main()` - 7 edges
8. `read_qc_summary()` - 6 edges
9. `main()` - 6 edges
10. `normalize_transcript()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `GALBA2` --conceptually_related_to--> `GALBA2 logo: broccoli Roman centurion character`  [INFERRED]
  README.md → img/logo.png
- `GALBA2` --conceptually_related_to--> `GALBA2 pipeline overview diagram`  [INFERRED]
  README.md → img/pipeline_overview.png
- `barrnap (rRNA predictor)` --references--> `GALBA2 pipeline overview diagram`  [EXTRACTED]
  README.md → img/pipeline_overview.png
- `tRNAscan-SE (tRNA predictor)` --references--> `GALBA2 pipeline overview diagram`  [EXTRACTED]
  README.md → img/pipeline_overview.png
- `Infernal + Rfam (ncRNA scanner)` --references--> `GALBA2 pipeline overview diagram`  [EXTRACTED]
  README.md → img/pipeline_overview.png

## Hyperedges (group relationships)
- **Core GALBA2 Pipeline Tools** — readme_miniprot, readme_miniprothint, readme_augustus, readme_agat, readme_diamond [EXTRACTED 1.00]
- **Optional ncRNA Annotation Tools** — readme_barrnap, readme_trnascan_se, readme_infernal [EXTRACTED 1.00]
- **Optional QC Tools** — readme_busco, readme_compleasm, readme_omark, readme_gffcompare [EXTRACTED 1.00]
- **Repeat Masking Tools (alternative)** — readme_repeatmodeler2, readme_repeatmasker, readme_red [EXTRACTED 1.00]
- **GALBA2 Input Configuration Files** — readme_samples_csv, readme_config_ini [EXTRACTED 1.00]
- **Host-side Python Dependencies** — requirements_snakemake, requirements_slurm_plugin, requirements_pandas [EXTRACTED 1.00]
- **Pipeline Overview Visualizations** — pipeline_overview_diagram, pipeline_overview_svg [EXTRACTED 1.00]
- **GTF Manipulation Scripts** — filter_gtf_by_txid_main, filter_gtf_by_diamond_main, normalize_gtf_main, gtf_sanity_check_main, get_longest_isoform_main, fix_in_frame_stop_main [INFERRED 0.90]
- **Hints File Generation** — generate_galba_hints_main, compleasm_to_hints_main, analyze_hint_support_main, gene_support_summary_main [INFERRED 0.90]
- **QC and Reporting** — generate_report_main, training_summary_main, completeness_plot_main, gene_set_statistics_main [INFERRED 0.90]
- **FANTASIA Functional Annotation** — fantasia_summary_main, fantasia_decorate_gff3_main [INFERRED 0.90]
- **ncRNA Prediction** — infernal_to_gff3_main [INFERRED 0.85]
- **Protein/Coding Sequence Extraction** — getannofasta_script, normalize_gtf_main, fix_in_frame_stop_main [INFERRED 0.90]
- **Repeat Masking Support** — parsetrf_main [INFERRED 0.85]

## Communities

### Community 0 - "Community 0"
Cohesion: 0.06
Nodes (51): Friends image: GALBA2 broccoli + hammer/bird character, GALBA2 logo: broccoli Roman centurion character, Deprecated: --crf/--keepCrf CRF training not ported, Deprecated: galba.pl --prg=gth (GenomeThreader) not ported, galba.pl --genome flag → samples.csv genome column, galba.pl --prot_seq flag → samples.csv protein_fasta column, Migration guide: galba.pl to GALBA2, GALBA2 pipeline overview diagram (+43 more)

### Community 1 - "Community 1"
Cohesion: 0.1
Nodes (31): collect_benchmarks(), deduplicate_bibtex(), deduplicate_citations(), detect_mode(), embed_image(), format_bbc_decisions(), format_compleasm_as_busco(), format_time() (+23 more)

### Community 2 - "Community 2"
Cohesion: 0.16
Nodes (17): analyze_hint_support(), check_cds_overlap(), check_intron_support(), load_hints(), load_transcripts_and_introns(), main(), parse_gff_line(), parse_gtf_line() (+9 more)

### Community 3 - "Community 3"
Cohesion: 0.17
Nodes (14): generate_plot(), main(), parse_busco_summary(), parse_compleasm_summary(), Parse BUSCO summary for genome and proteome scores., Parse compleasm summary.txt., Generate horizontal stacked bar chart., count_loci_in_gb() (+6 more)

### Community 4 - "Community 4"
Cohesion: 0.19
Nodes (15): get_cds_features(), get_sequence(), main(), normalize_transcript(), parse_attrs(), parse_gtf(), Update gene and transcript boundaries to match their features., Write normalized GTF, preserving gene/transcript/feature order. (+7 more)

### Community 5 - "Community 5"
Cohesion: 0.19
Nodes (11): check_exon_support(), check_support(), main(), parse_gtf_transcripts(), parse_hints(), Check how many intervals are supported by hints.      For introns: exact match (, Check how many exons have any overlapping exon/CDS hints.      Uses overlap (not, Parse hints file into interval trees by (chrom, feature_type, source_class). (+3 more)

### Community 6 - "Community 6"
Cohesion: 0.14
Nodes (10): check_tool_in_given_path(), create_log_file_name(), create_random_string(), create_tmp_dir(), find_tool(), Funtion that creates a random string added to the logfile name         and tmp d, Function that creates a log file with a random name, Function that creates a directory for temporary files with a random name (+2 more)

### Community 7 - "Community 7"
Cohesion: 0.27
Nodes (9): Enum, addToDict(), addToGc(), initCounts(), main(), parse(), parseCmd(), printStatistics() (+1 more)

### Community 8 - "Community 8"
Cohesion: 0.24
Nodes (11): compute_statistics(), generate_plots(), main(), parse_gtf(), parse_support_tsv(), Parse gene_support.tsv for evidence support visualization., Generate all publication-quality plots., Parse GTF file into gene/transcript/exon structure.      Returns:         genes: (+3 more)

### Community 9 - "Community 9"
Cohesion: 0.26
Nodes (11): classify_proteins(), _compile_categories(), main(), parse_results(), Map each protein to all categories its GO terms touch.      Returns a Counter of, Stream results.csv: collect counters AND per-row tuples for the TSV., Flat per-(transcript, GO term) TSV with human-readable GO names., Pie chart of broad functional categories.      Each protein is counted once per (+3 more)

### Community 10 - "Community 10"
Cohesion: 0.29
Nodes (9): _add_ontology_term(), decorate(), load_go_assignments(), _lookup_key(), main(), Two-pass: build transcript->gene map, then rewrite the GFF3.      Both `mRNA` an, Return the FANTASIA lookup key for a transcript-equivalent row.      Prefers `tr, Map transcript_id -> set of GO IDs above the score threshold. (+1 more)

### Community 11 - "Community 11"
Cohesion: 0.43
Nodes (6): check_gtf(), main(), Reads a GTF file and extracts CDS (Coding Sequence) information for each transcr, Checks the consistency and integrity of CDS (Coding Sequence) information within, read_gtf(), write_output()

### Community 12 - "Community 12"
Cohesion: 0.6
Nodes (5): extract_tx_ids_from_tsv(), main(), miniprot_to_hints(), read_and_filter_gff(), run_simple_process()

### Community 13 - "Community 13"
Cohesion: 0.47
Nodes (5): filter_gtf(), main(), Run DIAMOND to search for homologous proteins in the reference protein set, Read DIAMOND output, identify transcript names that have at least one hit,     r, run_diamond()

### Community 14 - "Community 14"
Cohesion: 0.6
Nodes (4): main(), parse_tblout(), Parse Infernal --fmt 2 tblout output.      Fields (fmt 2):     0: idx, 1: target, write_gff3()

### Community 15 - "Community 15"
Cohesion: 1.0
Nodes (0): 

### Community 16 - "Community 16"
Cohesion: 1.0
Nodes (2): cfg/galba.cfg, config: galba_cfg_path

### Community 17 - "Community 17"
Cohesion: 1.0
Nodes (0): 

### Community 18 - "Community 18"
Cohesion: 1.0
Nodes (0): 

### Community 19 - "Community 19"
Cohesion: 1.0
Nodes (1): config: samples_file

### Community 20 - "Community 20"
Cohesion: 1.0
Nodes (1): config: min_contig

### Community 21 - "Community 21"
Cohesion: 1.0
Nodes (1): config: slurm_args

### Community 22 - "Community 22"
Cohesion: 1.0
Nodes (1): config: masking_tool

### Community 23 - "Community 23"
Cohesion: 1.0
Nodes (1): config: skip_busco

### Community 24 - "Community 24"
Cohesion: 1.0
Nodes (1): config: run_omark

### Community 25 - "Community 25"
Cohesion: 1.0
Nodes (1): config: run_ncrna

### Community 26 - "Community 26"
Cohesion: 1.0
Nodes (1): config: run_fantasia

### Community 27 - "Community 27"
Cohesion: 1.0
Nodes (1): config: disable_diamond_filter

### Community 28 - "Community 28"
Cohesion: 1.0
Nodes (1): config: translation_table

### Community 29 - "Community 29"
Cohesion: 1.0
Nodes (1): container: galba-tools (katharinahoff/galba2-tools:latest)

### Community 30 - "Community 30"
Cohesion: 1.0
Nodes (1): container: augustus (quay.io/biocontainers/augustus:3.5.0)

### Community 31 - "Community 31"
Cohesion: 1.0
Nodes (1): container: AGAT (quay.io/biocontainers/agat:1.4.1)

### Community 32 - "Community 32"
Cohesion: 1.0
Nodes (1): container: Red (quay.io/biocontainers/red:2018.09.10)

### Community 33 - "Community 33"
Cohesion: 1.0
Nodes (1): container: dfam/tetools:latest (RepeatModeler+RepeatMasker)

### Community 34 - "Community 34"
Cohesion: 1.0
Nodes (1): file: output/{sample}/genome.fa

### Community 35 - "Community 35"
Cohesion: 1.0
Nodes (1): file: output/{sample}/genome_masked.fa

### Community 36 - "Community 36"
Cohesion: 1.0
Nodes (1): file: output/{sample}/preprocessing/genome.fa.masked

### Community 37 - "Community 37"
Cohesion: 1.0
Nodes (1): file: output/{sample}/preprocessing/proteins_merged.fa

### Community 38 - "Community 38"
Cohesion: 1.0
Nodes (1): file: output/{sample}/preprocessing/.headers_fixed

### Community 39 - "Community 39"
Cohesion: 1.0
Nodes (1): file: output/{sample}/augustus.hints.gtf

### Community 40 - "Community 40"
Cohesion: 1.0
Nodes (1): file: output/{sample}/bad_genes.lst

### Community 41 - "Community 41"
Cohesion: 1.0
Nodes (1): file: output/{sample}/augustus.hints.fixed.gtf

### Community 42 - "Community 42"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.filtered.gtf

### Community 43 - "Community 43"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.normalized.gtf

### Community 44 - "Community 44"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.gtf

### Community 45 - "Community 45"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.fixed.gtf

### Community 46 - "Community 46"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.gff3

### Community 47 - "Community 47"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.aa

### Community 48 - "Community 48"
Cohesion: 1.0
Nodes (1): file: output/{sample}/galba.codingseq

### Community 49 - "Community 49"
Cohesion: 1.0
Nodes (1): file: output/{sample}/bonafide.f.clean.gb

### Community 50 - "Community 50"
Cohesion: 1.0
Nodes (1): file: output/{sample}/miniprot/miniprot_trainingGenes.gtf

### Community 51 - "Community 51"
Cohesion: 1.0
Nodes (1): file: output/{sample}/bonafide.f.gb

### Community 52 - "Community 52"
Cohesion: 1.0
Nodes (1): file: output/{sample}/bonafide.f.gtf

### Community 53 - "Community 53"
Cohesion: 1.0
Nodes (1): file: output/{sample}/hintsfile.gff

### Community 54 - "Community 54"
Cohesion: 1.0
Nodes (1): file: output/{sample}/fantasia/results.csv

### Community 55 - "Community 55"
Cohesion: 1.0
Nodes (1): file: output/{sample}/fantasia/galba.go.gff3

### Community 56 - "Community 56"
Cohesion: 1.0
Nodes (1): file: output/{sample}/results/.done

### Community 57 - "Community 57"
Cohesion: 1.0
Nodes (1): tool: RepeatModeler2

### Community 58 - "Community 58"
Cohesion: 1.0
Nodes (1): tool: RepeatMasker

### Community 59 - "Community 59"
Cohesion: 1.0
Nodes (1): tool: TRF (Tandem Repeat Finder)

### Community 60 - "Community 60"
Cohesion: 1.0
Nodes (1): tool: Red (REpeat Detector)

### Community 61 - "Community 61"
Cohesion: 1.0
Nodes (1): tool: agat_convert_sp_gxf2gxf.pl

### Community 62 - "Community 62"
Cohesion: 1.0
Nodes (1): tool: DIAMOND

### Community 63 - "Community 63"
Cohesion: 1.0
Nodes (1): tool: aa2nonred.pl

### Community 64 - "Community 64"
Cohesion: 1.0
Nodes (1): tool: filterGenesIn.pl

### Community 65 - "Community 65"
Cohesion: 1.0
Nodes (1): tool: singularity (FANTASIA-Lite GPU container)

### Community 66 - "Community 66"
Cohesion: 1.0
Nodes (1): tool: compleasm.py

### Community 67 - "Community 67"
Cohesion: 1.0
Nodes (1): barrnap

### Community 68 - "Community 68"
Cohesion: 1.0
Nodes (1): tRNAscan-SE

### Community 69 - "Community 69"
Cohesion: 1.0
Nodes (1): Infernal/cmscan

### Community 70 - "Community 70"
Cohesion: 1.0
Nodes (1): AUGUSTUS

### Community 71 - "Community 71"
Cohesion: 1.0
Nodes (1): miniprot

### Community 72 - "Community 72"
Cohesion: 1.0
Nodes (1): minisplice

### Community 73 - "Community 73"
Cohesion: 1.0
Nodes (1): BUSCO

### Community 74 - "Community 74"
Cohesion: 1.0
Nodes (1): gffcompare

### Community 75 - "Community 75"
Cohesion: 1.0
Nodes (1): OMArk

## Knowledge Gaps
- **131 isolated node(s):** `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.`, `Read file contents, return empty string if not found.`, `Base64-encode an image for inline HTML embedding.      If download_name is provi` (+126 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 15`** (2 nodes): `config.ini`, `samples.csv`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 16`** (2 nodes): `cfg/galba.cfg`, `config: galba_cfg_path`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 17`** (1 nodes): `filter_gtf_by_txid.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 18`** (1 nodes): `getAnnoFastaFromJoingenes.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 19`** (1 nodes): `config: samples_file`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 20`** (1 nodes): `config: min_contig`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 21`** (1 nodes): `config: slurm_args`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 22`** (1 nodes): `config: masking_tool`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 23`** (1 nodes): `config: skip_busco`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 24`** (1 nodes): `config: run_omark`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 25`** (1 nodes): `config: run_ncrna`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 26`** (1 nodes): `config: run_fantasia`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 27`** (1 nodes): `config: disable_diamond_filter`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 28`** (1 nodes): `config: translation_table`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 29`** (1 nodes): `container: galba-tools (katharinahoff/galba2-tools:latest)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 30`** (1 nodes): `container: augustus (quay.io/biocontainers/augustus:3.5.0)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 31`** (1 nodes): `container: AGAT (quay.io/biocontainers/agat:1.4.1)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 32`** (1 nodes): `container: Red (quay.io/biocontainers/red:2018.09.10)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 33`** (1 nodes): `container: dfam/tetools:latest (RepeatModeler+RepeatMasker)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 34`** (1 nodes): `file: output/{sample}/genome.fa`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 35`** (1 nodes): `file: output/{sample}/genome_masked.fa`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 36`** (1 nodes): `file: output/{sample}/preprocessing/genome.fa.masked`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 37`** (1 nodes): `file: output/{sample}/preprocessing/proteins_merged.fa`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 38`** (1 nodes): `file: output/{sample}/preprocessing/.headers_fixed`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 39`** (1 nodes): `file: output/{sample}/augustus.hints.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 40`** (1 nodes): `file: output/{sample}/bad_genes.lst`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 41`** (1 nodes): `file: output/{sample}/augustus.hints.fixed.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 42`** (1 nodes): `file: output/{sample}/galba.filtered.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 43`** (1 nodes): `file: output/{sample}/galba.normalized.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 44`** (1 nodes): `file: output/{sample}/galba.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 45`** (1 nodes): `file: output/{sample}/galba.fixed.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 46`** (1 nodes): `file: output/{sample}/galba.gff3`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 47`** (1 nodes): `file: output/{sample}/galba.aa`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 48`** (1 nodes): `file: output/{sample}/galba.codingseq`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 49`** (1 nodes): `file: output/{sample}/bonafide.f.clean.gb`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 50`** (1 nodes): `file: output/{sample}/miniprot/miniprot_trainingGenes.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 51`** (1 nodes): `file: output/{sample}/bonafide.f.gb`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 52`** (1 nodes): `file: output/{sample}/bonafide.f.gtf`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 53`** (1 nodes): `file: output/{sample}/hintsfile.gff`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 54`** (1 nodes): `file: output/{sample}/fantasia/results.csv`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 55`** (1 nodes): `file: output/{sample}/fantasia/galba.go.gff3`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 56`** (1 nodes): `file: output/{sample}/results/.done`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 57`** (1 nodes): `tool: RepeatModeler2`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 58`** (1 nodes): `tool: RepeatMasker`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 59`** (1 nodes): `tool: TRF (Tandem Repeat Finder)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 60`** (1 nodes): `tool: Red (REpeat Detector)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 61`** (1 nodes): `tool: agat_convert_sp_gxf2gxf.pl`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 62`** (1 nodes): `tool: DIAMOND`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 63`** (1 nodes): `tool: aa2nonred.pl`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 64`** (1 nodes): `tool: filterGenesIn.pl`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 65`** (1 nodes): `tool: singularity (FANTASIA-Lite GPU container)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 66`** (1 nodes): `tool: compleasm.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 67`** (1 nodes): `barrnap`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 68`** (1 nodes): `tRNAscan-SE`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 69`** (1 nodes): `Infernal/cmscan`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 70`** (1 nodes): `AUGUSTUS`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 71`** (1 nodes): `miniprot`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 72`** (1 nodes): `minisplice`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 73`** (1 nodes): `BUSCO`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 74`** (1 nodes): `gffcompare`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 75`** (1 nodes): `OMArk`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `main()` connect `Community 2` to `Community 5`?**
  _High betweenness centrality (0.008) - this node is a cross-community bridge._
- **Why does `main()` connect `Community 5` to `Community 2`?**
  _High betweenness centrality (0.007) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `GALBA2` (e.g. with `GALBA2 logo: broccoli Roman centurion character` and `GALBA2 pipeline overview diagram`) actually correct?**
  _`GALBA2` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `main()` (e.g. with `main()` and `main()`) actually correct?**
  _`main()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.` to the rest of the system?**
  _131 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.06 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.1 - nodes in this community are weakly interconnected._