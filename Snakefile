"""
GALBA2 Snakemake Workflow

Protein-only gene prediction pipeline using:
- Miniprot for protein-to-genome alignment
- AUGUSTUS for gene prediction with protein-derived hints

Authors: Katharina J. Hoff
Version: 0.1.0-beta
"""

__author__ = "Katharina J. Hoff"
__version__ = "0.1.0-beta"

import pandas as pd
import configparser
import os

# ============================================================================
# Configuration
# ============================================================================

# Read config.ini
config_ini_path = os.environ.get('GALBA2_CONFIG', 'config.ini')
# Allow trailing # / ; comments on value lines (ConfigParser default is None).
config_parser = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
if not os.path.isfile(config_ini_path):
    raise FileNotFoundError(
        f"config.ini not found at '{config_ini_path}'. "
        "Set GALBA2_CONFIG to a valid path or create config.ini in the "
        "working directory. See the README for an example configuration."
    )
config_parser.read(config_ini_path)

# Overlay environment-variable overrides
for section in ('PARAMS', 'SLURM_ARGS', 'fantasia'):
    if not config_parser.has_section(section):
        config_parser.add_section(section)

_env_overrides = {
    # SLURM_ARGS
    'GALBA2_CPUS_PER_TASK':                  ('SLURM_ARGS', 'cpus_per_task'),
    'GALBA2_MEM_OF_NODE':                    ('SLURM_ARGS', 'mem_of_node'),
    'GALBA2_MAX_RUNTIME':                    ('SLURM_ARGS', 'max_runtime'),
    # PARAMS
    'GALBA2_MIN_CONTIG':                     ('PARAMS', 'min_contig'),
    'GALBA2_SKIP_OPTIMIZE_AUGUSTUS':          ('PARAMS', 'skip_optimize_augustus'),
    'GALBA2_USE_DEV_SHM':                    ('PARAMS', 'use_dev_shm'),
    'GALBA2_SKIP_BUSCO':                     ('PARAMS', 'skip_busco'),
    'GALBA2_RUN_OMARK':                      ('PARAMS', 'run_omark'),
    'GALBA2_NO_CLEANUP':                     ('PARAMS', 'no_cleanup'),
    'GALBA2_TRANSLATION_TABLE':              ('PARAMS', 'translation_table'),
    'GALBA2_GC_DONOR':                       ('PARAMS', 'gc_donor'),
    'GALBA2_ALLOW_HINTED_SPLICESITES':       ('PARAMS', 'allow_hinted_splicesites'),
    'GALBA2_DISABLE_DIAMOND_FILTER':         ('PARAMS', 'disable_diamond_filter'),
    'GALBA2_AUGUSTUS_CHUNKSIZE':              ('PARAMS', 'augustus_chunksize'),
    'GALBA2_AUGUSTUS_OVERLAP':                ('PARAMS', 'augustus_overlap'),
    # FANTASIA-Lite (optional, GPU-only functional annotation)
    'GALBA2_RUN_FANTASIA':                   ('fantasia', 'enable'),
    'GALBA2_FANTASIA_SIF':                   ('fantasia', 'sif'),
    'GALBA2_FANTASIA_HF_CACHE':              ('fantasia', 'hf_cache_dir'),
    'GALBA2_FANTASIA_PARTITION':             ('fantasia', 'partition'),
    'GALBA2_FANTASIA_GPUS':                  ('fantasia', 'gpus'),
    'GALBA2_FANTASIA_MEM_MB':                ('fantasia', 'mem_mb'),
    'GALBA2_FANTASIA_CPUS':                  ('fantasia', 'cpus_per_task'),
    'GALBA2_FANTASIA_MAX_RUNTIME':           ('fantasia', 'max_runtime'),
    'GALBA2_FANTASIA_MIN_SCORE':             ('fantasia', 'min_score'),
    'GALBA2_FANTASIA_ADDITIONAL_PARAMS':     ('fantasia', 'additional_params'),
}
for _env_name, (_section, _key) in _env_overrides.items():
    if _env_name in os.environ:
        config_parser.set(_section, _key, os.environ[_env_name])

augustus_config_path = os.path.abspath(config_parser['paths'].get('augustus_config_path', 'augustus_config'))

# Pass config to Snakemake config dict
config['samples_file'] = config_parser['paths'].get('samples_file', 'samples.csv')
config['min_contig'] = config_parser.getint('PARAMS', 'min_contig', fallback=10000)
config['galba_tools_image'] = config_parser.get('containers', 'galba_tools_image',
                                                fallback='docker://katharinahoff/galba-notebook:latest')
config['augustus_image'] = config_parser.get('containers', 'augustus_image',
                                            fallback='docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9716f88_9')
# Backwards compat alias
config['galba_image'] = config['galba_tools_image']

config['slurm_args'] = {
    'cpus_per_task': config_parser.getint('SLURM_ARGS', 'cpus_per_task', fallback=1),
    'mem_of_node': config_parser.getint('SLURM_ARGS', 'mem_of_node', fallback=16000),
    'max_runtime': config_parser.getint('SLURM_ARGS', 'max_runtime', fallback=60)
}

config['skip_optimize_augustus'] = config_parser.getboolean(
    'PARAMS', 'skip_optimize_augustus', fallback=False
)

config['use_dev_shm'] = config_parser.getboolean(
    'PARAMS', 'use_dev_shm', fallback=False
)

config['skip_busco'] = config_parser.getboolean(
    'PARAMS', 'skip_busco', fallback=False
)

config['run_omark'] = config_parser.getboolean(
    'PARAMS', 'run_omark', fallback=False
)

config['no_cleanup'] = config_parser.getboolean(
    'PARAMS', 'no_cleanup', fallback=False
)

config['translation_table'] = config_parser.getint(
    'PARAMS', 'translation_table', fallback=1
)

config['gc_donor'] = config_parser.getfloat(
    'PARAMS', 'gc_donor', fallback=0.001
)

config['allow_hinted_splicesites'] = config_parser.get(
    'PARAMS', 'allow_hinted_splicesites', fallback='gcag,atac'
)

config['disable_diamond_filter'] = config_parser.getboolean(
    'PARAMS', 'disable_diamond_filter', fallback=False
)

config['augustus_chunksize'] = config_parser.getint(
    'PARAMS', 'augustus_chunksize', fallback=3000000
)
config['augustus_overlap'] = config_parser.getint(
    'PARAMS', 'augustus_overlap', fallback=500000
)

# FANTASIA-Lite functional annotation (optional, GPU-only).
config['run_fantasia'] = config_parser.getboolean(
    'fantasia', 'enable', fallback=False
)
config['fantasia'] = {
    'sif':              config_parser.get('fantasia', 'sif', fallback=''),
    'hf_cache_dir':     config_parser.get('fantasia', 'hf_cache_dir', fallback=''),
    'additional_params': config_parser.get('fantasia', 'additional_params', fallback=''),
    'min_score':        config_parser.get('fantasia', 'min_score', fallback='0.5'),
    'partition':        config_parser.get('fantasia', 'partition', fallback=''),
    'gpus':             config_parser.get('fantasia', 'gpus', fallback='1'),
    'mem_mb':           config_parser.get('fantasia', 'mem_mb',
                                          fallback=str(config['slurm_args']['mem_of_node'])),
    'cpus_per_task':    config_parser.get('fantasia', 'cpus_per_task',
                                          fallback=str(config['slurm_args']['cpus_per_task'])),
    'max_runtime':      config_parser.get('fantasia', 'max_runtime',
                                          fallback=str(config['slurm_args']['max_runtime'])),
}

config['username'] = os.environ.get("USER", "unknown")
script_dir = os.path.join(os.path.dirname(workflow.main_snakefile), "scripts")

# Shared data paths — auto-discover from this repo, sibling BRAKER4/BRAKER-as-snakemake,
# or use explicit paths from config.ini.
def _find_shared_data():
    """Find shared_data directory: this repo, then sibling repos, then create."""
    repo_dir = os.path.dirname(workflow.main_snakefile)
    # 1. This repo
    local = os.path.join(repo_dir, 'shared_data')
    if os.path.isdir(local):
        return local
    # 2. Sibling repos
    parent = os.path.dirname(repo_dir)
    for sibling in ('BRAKER4', 'BRAKER-as-snakemake'):
        sib_shared = os.path.join(parent, sibling, 'shared_data')
        if os.path.isdir(sib_shared):
            return sib_shared
    # 3. Fallback
    return local

_shared_data = _find_shared_data()
config['busco_download_path'] = os.path.abspath(
    config_parser.get('paths', 'busco_download_path',
                      fallback=os.path.join(_shared_data, 'busco_downloads'))
)
config['compleasm_download_path'] = os.path.abspath(
    config_parser.get('paths', 'compleasm_download_path',
                      fallback=os.path.join(_shared_data, 'compleasm_downloads'))
)

# OMAmer database path (for optional OMArk QC)
def _find_omamer_db():
    """Find LUCA.h5: config, shared_data, sibling repos."""
    explicit = config_parser.get('PARAMS', 'omamer_db', fallback='')
    if explicit and os.path.isfile(explicit):
        return explicit
    # Check shared_data
    for candidate in [
        os.path.join(_shared_data, 'LUCA.h5'),
        os.path.join(os.path.dirname(workflow.main_snakefile), 'test_data', 'LUCA.h5'),
    ]:
        if os.path.isfile(candidate):
            return candidate
    # Check sibling repos
    parent = os.path.dirname(os.path.dirname(workflow.main_snakefile))
    for sibling in ('BRAKER4', 'BRAKER-as-snakemake'):
        for subdir in ('shared_data', 'test_data'):
            candidate = os.path.join(parent, sibling, subdir, 'LUCA.h5')
            if os.path.isfile(candidate):
                return candidate
    return os.path.join(_shared_data, 'LUCA.h5')

config['omamer_db'] = _find_omamer_db()

# Path to the GALBA extrinsic config file bundled with the pipeline
config['galba_cfg_path'] = os.path.join(
    os.path.dirname(workflow.main_snakefile), "cfg", "galba.cfg"
)

config['pipeline_version'] = __version__

# ============================================================================
# Include Common Functions and Sample Parsing
# ============================================================================

include: "rules/common.smk"

# ============================================================================
# Detect Data Types
# ============================================================================

GLOBAL_DATA_TYPES = {
    'has_proteins': False,
    'has_reference_gtf': False,
    'needs_masking': False,
}

for sample in SAMPLES:
    types = detect_data_types(sample)
    if types['has_proteins']:
        GLOBAL_DATA_TYPES['has_proteins'] = True
    if types.get('has_reference_gtf'):
        GLOBAL_DATA_TYPES['has_reference_gtf'] = True
    if types.get('needs_masking'):
        GLOBAL_DATA_TYPES['needs_masking'] = True

print("\nGALBA2 pipeline rules:")
print("  ✓ Genome preparation (header cleaning)")
if GLOBAL_DATA_TYPES['needs_masking']:
    print("  ✓ RepeatModeler2 + RepeatMasker (genome masking)")
print("  ✓ Miniprot protein-to-genome alignment")
print("  ✓ Miniprothint training gene and hint generation")
print("  ✓ AUGUSTUS training on protein-derived genes")
print("  ✓ AUGUSTUS prediction with protein hints")
if not config.get('disable_diamond_filter', False):
    print("  ✓ DIAMOND filtering of predictions against input proteins")
if config.get('run_fantasia', False):
    print("  ✓ FANTASIA-Lite functional annotation (GPU)")
if GLOBAL_DATA_TYPES['has_reference_gtf']:
    print("  ✓ gffcompare evaluation against reference annotation")
print()

# ============================================================================
# Include Rules
# ============================================================================

# Preprocessing
include: "rules/preprocessing/prepare_genome.smk"
include: "rules/preprocessing/merge_proteins.smk"
if GLOBAL_DATA_TYPES['needs_masking']:
    include: "rules/preprocessing/run_masking.smk"

# Miniprot alignment and hint generation
include: "rules/miniprot/run_miniprot.smk"
include: "rules/miniprot/generate_hints.smk"

# AUGUSTUS training
include: "rules/augustus_training/copy_augustus_config.smk"
include: "rules/augustus_training/create_new_species.smk"
include: "rules/augustus_training/compute_flanking_region.smk"
include: "rules/augustus_training/convert_to_genbank.smk"
include: "rules/augustus_training/filter_genes_etraining.smk"
include: "rules/postprocessing/redundancy_removal.smk"
include: "rules/augustus_training/downsample_training_genes.smk"
include: "rules/augustus_training/split_training_set.smk"
include: "rules/augustus_training/train_augustus.smk"
include: "rules/augustus_training/optimize_augustus.smk"

# AUGUSTUS prediction
include: "rules/augustus_predict/run_augustus_hints.smk"

# Postprocessing
include: "rules/postprocessing/fix_in_frame_stop_codons.smk"
include: "rules/postprocessing/filter_stop_codons.smk"
include: "rules/postprocessing/normalize_cds.smk"
include: "rules/postprocessing/postprocess_augustus.smk"
include: "rules/postprocessing/convert_to_gff3.smk"

# Quality control
include: "rules/quality_control/run_compleasm.smk"
include: "rules/quality_control/gene_support.smk"
if not config.get('skip_busco', False):
    include: "rules/quality_control/run_busco.smk"
if config.get('run_omark', False):
    include: "rules/quality_control/run_omark.smk"
if GLOBAL_DATA_TYPES['has_reference_gtf']:
    include: "rules/quality_control/run_gffcompare.smk"

# FANTASIA-Lite functional annotation (optional, GPU-only)
if config.get('run_fantasia', False):
    include: "rules/postprocessing/run_fantasia.smk"

# Results collection
include: "rules/postprocessing/collect_results.smk"

# ============================================================================
# Target Rule
# ============================================================================

def get_final_outputs():
    """Determine final outputs based on sample configurations."""
    outputs = []
    for sample in SAMPLES:
        outputs.append(f"output/{sample}/results/.done")
    return outputs

rule all:
    input:
        get_final_outputs()
    default_target: True

# ============================================================================
# Wildcard Constraints
# ============================================================================

wildcard_constraints:
    sample="|".join(SAMPLES)

# ============================================================================
# Workflow Summary
# ============================================================================

print("\nTarget outputs:")
for output in get_final_outputs():
    print(f"  → {output}")
print(f"\n{'='*70}\n")
