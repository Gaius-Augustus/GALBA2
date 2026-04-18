# Graph Report - .  (2026-04-18)

## Corpus Check
- Large corpus: 1507 files · ~537,237 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder, or use --no-semantic to run AST-only.

## Summary
- 386 nodes · 560 edges · 24 communities detected
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 31 edges (avg confidence: 0.85)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Pipeline Orchestration & Data Flow|Pipeline Orchestration & Data Flow]]
- [[_COMMUNITY_GALBA2 Docs & Project Identity|GALBA2 Docs & Project Identity]]
- [[_COMMUNITY_HTML Report Generator|HTML Report Generator]]
- [[_COMMUNITY_AUGUSTUS Training & QC|AUGUSTUS Training & QC]]
- [[_COMMUNITY_Hint Support Analysis|Hint Support Analysis]]
- [[_COMMUNITY_FANTASIA GO Annotation|FANTASIA GO Annotation]]
- [[_COMMUNITY_Repeat Masking|Repeat Masking]]
- [[_COMMUNITY_BUSCO Completeness Plots|BUSCO Completeness Plots]]
- [[_COMMUNITY_GTF Normalization|GTF Normalization]]
- [[_COMMUNITY_Stop Codon Fix|Stop Codon Fix]]
- [[_COMMUNITY_GTF Filtering & Isoform Selection|GTF Filtering & Isoform Selection]]
- [[_COMMUNITY_Gene Set Statistics|Gene Set Statistics]]
- [[_COMMUNITY_TRF Repeat Parsing|TRF Repeat Parsing]]
- [[_COMMUNITY_ncRNA Annotation Rules|ncRNA Annotation Rules]]
- [[_COMMUNITY_compleasm to Hints|compleasm to Hints]]
- [[_COMMUNITY_BUSCO QC Rules|BUSCO QC Rules]]
- [[_COMMUNITY_DIAMOND Reference Filter|DIAMOND Reference Filter]]
- [[_COMMUNITY_Infernal ncRNA to GFF3|Infernal ncRNA to GFF3]]
- [[_COMMUNITY_Masked Genome Prep|Masked Genome Prep]]
- [[_COMMUNITY_GTF Transcript Filter|GTF Transcript Filter]]
- [[_COMMUNITY_Joingenes FASTA Extractor|Joingenes FASTA Extractor]]
- [[_COMMUNITY_Samples Config|Samples Config]]
- [[_COMMUNITY_Min Contig Config|Min Contig Config]]
- [[_COMMUNITY_Empty Hints Creation|Empty Hints Creation]]

## God Nodes (most connected - your core abstractions)
1. `GALBA2` - 43 edges
2. `rule: collect_results` - 19 edges
3. `main()` - 17 edges
4. `GALBA2 pipeline overview diagram` - 15 edges
5. `rule: sanity_check_augustus` - 13 edges
6. `rule: run_masking` - 12 edges
7. `rule: fix_in_frame_stop_codons` - 12 edges
8. `rule: redundancy_removal` - 11 edges
9. `main()` - 9 edges
10. `main()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `GALBA2 logo: broccoli Roman centurion character` --conceptually_related_to--> `GALBA2`  [INFERRED]
  img/logo.png → README.md
- `GALBA2 pipeline overview diagram` --conceptually_related_to--> `GALBA2`  [INFERRED]
  img/pipeline_overview.png → README.md
- `GALBA2 pipeline overview diagram` --references--> `barrnap (rRNA predictor)`  [EXTRACTED]
  img/pipeline_overview.png → README.md
- `GALBA2 pipeline overview diagram` --references--> `tRNAscan-SE (tRNA predictor)`  [EXTRACTED]
  img/pipeline_overview.png → README.md
- `GALBA2 pipeline overview diagram` --references--> `Infernal + Rfam (ncRNA scanner)`  [EXTRACTED]
  img/pipeline_overview.png → README.md

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

### Community 0 - "Pipeline Orchestration & Data Flow"
Cohesion: 0.05
Nodes (66): cfg/galba.cfg, rules/common.smk, config: disable_diamond_filter, config: galba_cfg_path, config: run_fantasia, config: run_ncrna, config: run_omark, config: skip_busco (+58 more)

### Community 1 - "GALBA2 Docs & Project Identity"
Cohesion: 0.06
Nodes (54): Friends image: GALBA2 broccoli + hammer/bird character, GALBA2 logo: broccoli Roman centurion character, Deprecated: --crf/--keepCrf CRF training not ported, Deprecated: galba.pl --prg=gth (GenomeThreader) not ported, galba.pl --genome flag → samples.csv genome column, galba.pl --prot_seq flag → samples.csv protein_fasta column, Migration guide: galba.pl to GALBA2, GALBA2 pipeline overview diagram (+46 more)

### Community 2 - "HTML Report Generator"
Cohesion: 0.1
Nodes (31): collect_benchmarks(), deduplicate_bibtex(), deduplicate_citations(), detect_mode(), embed_image(), format_bbc_decisions(), format_compleasm_as_busco(), format_time() (+23 more)

### Community 3 - "AUGUSTUS Training & QC"
Cohesion: 0.09
Nodes (29): cleanup_tmp_opt, compute_flanking_region, convert_to_genbank, copy_augustus_config, create_new_species, downsample_training_genes, filter_genes_etraining, gene_support_summary (+21 more)

### Community 4 - "Hint Support Analysis"
Cohesion: 0.1
Nodes (27): analyze_hint_support(), check_cds_overlap(), check_intron_support(), load_hints(), load_transcripts_and_introns(), main(), parse_gff_line(), parse_gtf_line() (+19 more)

### Community 5 - "FANTASIA GO Annotation"
Cohesion: 0.12
Nodes (22): _add_ontology_term(), decorate(), load_go_assignments(), _lookup_key(), main(), Two-pass: build transcript->gene map, then rewrite the GFF3.      Both `mRNA` an, Return the FANTASIA lookup key for a transcript-equivalent row.      Prefers `tr, Map transcript_id -> set of GO IDs above the score threshold. (+14 more)

### Community 6 - "Repeat Masking"
Cohesion: 0.16
Nodes (17): config: masking_tool, config: slurm_args, container: Red (quay.io/biocontainers/red:2018.09.10), container: dfam/tetools:latest (RepeatModeler+RepeatMasker), file: output/{sample}/genome.fa, file: output/{sample}/preprocessing/.headers_fixed, file: output/{sample}/preprocessing/genome.fa.masked, rule: prepare_genome (+9 more)

### Community 7 - "BUSCO Completeness Plots"
Cohesion: 0.17
Nodes (14): generate_plot(), main(), parse_busco_summary(), parse_compleasm_summary(), Parse BUSCO summary for genome and proteome scores., Parse compleasm summary.txt., Generate horizontal stacked bar chart., count_loci_in_gb() (+6 more)

### Community 8 - "GTF Normalization"
Cohesion: 0.19
Nodes (15): get_cds_features(), get_sequence(), main(), normalize_transcript(), parse_attrs(), parse_gtf(), Update gene and transcript boundaries to match their features., Write normalized GTF, preserving gene/transcript/feature order. (+7 more)

### Community 9 - "Stop Codon Fix"
Cohesion: 0.14
Nodes (10): check_tool_in_given_path(), create_log_file_name(), create_random_string(), create_tmp_dir(), find_tool(), Funtion that creates a random string added to the logfile name         and tmp d, Function that creates a log file with a random name, Function that creates a directory for temporary files with a random name (+2 more)

### Community 10 - "GTF Filtering & Isoform Selection"
Cohesion: 0.21
Nodes (9): filter_gtf_by_diamond_against_ref: main(), filter_gtf_by_txid: main (script-level), main(), check_gtf(), main(), Reads a GTF file and extracts CDS (Coding Sequence) information for each transcr, Checks the consistency and integrity of CDS (Coding Sequence) information within, read_gtf() (+1 more)

### Community 11 - "Gene Set Statistics"
Cohesion: 0.24
Nodes (11): compute_statistics(), generate_plots(), main(), parse_gtf(), parse_support_tsv(), Parse gene_support.tsv for evidence support visualization., Generate all publication-quality plots., Parse GTF file into gene/transcript/exon structure.      Returns:         genes: (+3 more)

### Community 12 - "TRF Repeat Parsing"
Cohesion: 0.27
Nodes (9): Enum, addToDict(), addToGc(), initCounts(), main(), parse(), parseCmd(), printStatistics() (+1 more)

### Community 13 - "ncRNA Annotation Rules"
Cohesion: 0.2
Nodes (9): merge_ncrna_into_gff3, merge_rrna_into_gff3, run_barrnap, convert_infernal_to_gff3, run_cmscan, run_trnascan, barrnap, Infernal/cmscan (+1 more)

### Community 14 - "compleasm to Hints"
Cohesion: 0.39
Nodes (6): extract_tx_ids_from_tsv(), main(), miniprot_to_hints(), read_and_filter_gff(), run_simple_process(), parseTrfOutput: main()

### Community 15 - "BUSCO QC Rules"
Cohesion: 0.33
Nodes (5): busco_genome, busco_proteins, busco_summary, get_longest_isoform, BUSCO

### Community 16 - "DIAMOND Reference Filter"
Cohesion: 0.47
Nodes (5): filter_gtf(), main(), Run DIAMOND to search for homologous proteins in the reference protein set, Read DIAMOND output, identify transcript names that have at least one hit,     r, run_diamond()

### Community 17 - "Infernal ncRNA to GFF3"
Cohesion: 0.6
Nodes (4): main(), parse_tblout(), Parse Infernal --fmt 2 tblout output.      Fields (fmt 2):     0: idx, 1: target, write_gff3()

### Community 18 - "Masked Genome Prep"
Cohesion: 1.0
Nodes (2): file: output/{sample}/genome_masked.fa, rule: prepare_masked_genome

### Community 19 - "GTF Transcript Filter"
Cohesion: 1.0
Nodes (0): 

### Community 20 - "Joingenes FASTA Extractor"
Cohesion: 1.0
Nodes (0): 

### Community 21 - "Samples Config"
Cohesion: 1.0
Nodes (1): config: samples_file

### Community 22 - "Min Contig Config"
Cohesion: 1.0
Nodes (1): config: min_contig

### Community 23 - "Empty Hints Creation"
Cohesion: 1.0
Nodes (1): rule: create_empty_hints

## Knowledge Gaps
- **123 isolated node(s):** `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.`, `Read file contents, return empty string if not found.`, `Base64-encode an image for inline HTML embedding.      If download_name is provi` (+118 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Masked Genome Prep`** (2 nodes): `file: output/{sample}/genome_masked.fa`, `rule: prepare_masked_genome`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `GTF Transcript Filter`** (1 nodes): `filter_gtf_by_txid.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Joingenes FASTA Extractor`** (1 nodes): `getAnnoFastaFromJoingenes.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Samples Config`** (1 nodes): `config: samples_file`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Min Contig Config`** (1 nodes): `config: min_contig`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Empty Hints Creation`** (1 nodes): `rule: create_empty_hints`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `main()` connect `HTML Report Generator` to `Infernal ncRNA to GFF3`, `Gene Set Statistics`, `FANTASIA GO Annotation`, `BUSCO Completeness Plots`?**
  _High betweenness centrality (0.097) - this node is a cross-community bridge._
- **Why does `main()` connect `Gene Set Statistics` to `HTML Report Generator`, `Hint Support Analysis`, `FANTASIA GO Annotation`?**
  _High betweenness centrality (0.089) - this node is a cross-community bridge._
- **Why does `main()` connect `Hint Support Analysis` to `GTF Filtering & Isoform Selection`, `Gene Set Statistics`?**
  _High betweenness centrality (0.075) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `GALBA2` (e.g. with `GALBA2 logo: broccoli Roman centurion character` and `GALBA2 pipeline overview diagram`) actually correct?**
  _`GALBA2` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Parse BUSCO summary for genome and proteome scores.`, `Parse compleasm summary.txt.`, `Generate horizontal stacked bar chart.` to the rest of the system?**
  _123 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Pipeline Orchestration & Data Flow` be split into smaller, more focused modules?**
  _Cohesion score 0.05 - nodes in this community are weakly interconnected._
- **Should `GALBA2 Docs & Project Identity` be split into smaller, more focused modules?**
  _Cohesion score 0.06 - nodes in this community are weakly interconnected._