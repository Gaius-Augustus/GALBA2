"""
GALBA2 Snakemake Workflow

Protein-only gene prediction pipeline using:
- Miniprot for protein-to-genome alignment
- AUGUSTUS for gene prediction with protein-derived hints

Authors: Katharina J. Hoff
Version: 0.1.0-beta
"""

__author__ = "Katharina J. Hoff"
__version__ = "0.4.0-beta"

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
for section in ('PARAMS', 'SLURM_ARGS', 'fantasia', 'OMARK'):
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
    'GALBA2_OMAMER_DB':                      ('OMARK', 'omamer_db'),
    'GALBA2_NO_CLEANUP':                     ('PARAMS', 'no_cleanup'),
    'GALBA2_TRANSLATION_TABLE':              ('PARAMS', 'translation_table'),
    'GALBA2_GC_DONOR':                       ('PARAMS', 'gc_donor'),
    'GALBA2_ALLOW_HINTED_SPLICESITES':       ('PARAMS', 'allow_hinted_splicesites'),
    'GALBA2_DISABLE_DIAMOND_FILTER':         ('PARAMS', 'disable_diamond_filter'),
    'GALBA2_AUGUSTUS_CHUNKSIZE':              ('PARAMS', 'augustus_chunksize'),
    'GALBA2_AUGUSTUS_OVERLAP':                ('PARAMS', 'augustus_overlap'),
    'GALBA2_MASKING_TOOL':                   ('PARAMS', 'masking_tool'),
    'GALBA2_USE_MINISPLICE':                 ('PARAMS', 'use_minisplice'),
    'GALBA2_MINISPLICE_MODEL':               ('PARAMS', 'minisplice_model'),
    'GALBA2_RUN_NCRNA':                      ('PARAMS', 'run_ncrna'),
    # FANTASIA-Lite (optional, GPU-only functional annotation)
    'GALBA2_RUN_FANTASIA':                   ('fantasia', 'enable'),
    'GALBA2_FANTASIA_SIF':                   ('fantasia', 'sif'),
    'GALBA2_FANTASIA_HF_CACHE':              ('fantasia', 'hf_cache_dir'),
    'GALBA2_FANTASIA_LOOKUP_DIR':            ('fantasia', 'lookup_dir'),
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
_container_defaults = {
    'galba_tools_image': 'docker://katharinahoff/galba-notebook:latest',
    'augustus_image':    'docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9716f88_9',
    'gffcompare_image':  'docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1',
    'red_image':         'docker://quay.io/biocontainers/red:2018.09.10--h9948957_3',
    'barrnap_image':     'docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4',
    'agat_image':        'docker://quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0',
    'busco_image':       'docker://ezlabgva/busco:v6.0.0_cv1',
    'omark_image':       'docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0',
    'tetools_image':     'docker://dfam/tetools:latest',
    'trnascan_image':    'docker://quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_0',
    'infernal_image':    'docker://quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_2',
}
for _img_key, _img_default in _container_defaults.items():
    config[_img_key] = config_parser.get('containers', _img_key, fallback=_img_default)
# Backwards compat alias
config['galba_image'] = config['galba_tools_image']

# When [SLURM_ARGS] cpus_per_task is missing from config.ini (typical for
# local runs), fall back to workflow.cores so that `snakemake --cores N`
# controls per-rule parallelism. Without this fallback, multithreaded rules
# would run single-threaded regardless of --cores.
# Same fix as BRAKER4 issue #10.
config['slurm_args'] = {
    'cpus_per_task': config_parser.getint('SLURM_ARGS', 'cpus_per_task',
                                          fallback=workflow.cores or 1),
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

config['masking_tool'] = config_parser.get(
    'PARAMS', 'masking_tool', fallback='repeatmasker'
)

config['use_minisplice'] = config_parser.getboolean(
    'PARAMS', 'use_minisplice', fallback=False
)
config['minisplice_model'] = config_parser.get(
    'PARAMS', 'minisplice_model', fallback=''
)

config['run_ncrna'] = config_parser.getboolean(
    'PARAMS', 'run_ncrna', fallback=False
)

# FANTASIA-Lite functional annotation (optional, GPU-only).
config['run_fantasia'] = config_parser.getboolean(
    'fantasia', 'enable', fallback=False
)
config['fantasia'] = {
    'sif':              config_parser.get('fantasia', 'sif', fallback=''),
    'hf_cache_dir':     config_parser.get('fantasia', 'hf_cache_dir', fallback=''),
    'lookup_dir':       config_parser.get('fantasia', 'lookup_dir', fallback=''),
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
if config['run_fantasia']:
    if not config['fantasia']['sif']:
        raise ValueError(
            "fantasia.enable=1 but fantasia.sif is empty. Set the path to the "
            "FANTASIA-Lite Singularity image (pre-pulled from "
            "docker://katharinahoff/fantasia_for_brain:lite.v1.0.0) in config.ini "
            "[fantasia] sif=, or via GALBA2_FANTASIA_SIF."
        )
    if not os.path.isfile(config['fantasia']['sif']):
        raise FileNotFoundError(
            f"fantasia.enable=1 but the configured SIF does not exist on disk: "
            f"{config['fantasia']['sif']}\n"
            "GALBA2 will not start a run that would only fail later inside the "
            "FANTASIA rule. Pre-pull the image with:\n"
            "    singularity pull <path> docker://katharinahoff/fantasia_for_brain:lite.v1.0.0"
        )
    if not config['fantasia']['hf_cache_dir']:
        raise ValueError(
            "fantasia.enable=1 but fantasia.hf_cache_dir is empty. Pre-cache "
            "the Rostlab/prot_t5_xl_uniref50 HuggingFace model and set "
            "fantasia.hf_cache_dir to the cache directory in config.ini, or via "
            "GALBA2_FANTASIA_HF_CACHE."
        )
    if not os.path.isdir(config['fantasia']['hf_cache_dir']):
        raise FileNotFoundError(
            f"fantasia.enable=1 but the configured HuggingFace cache directory "
            f"does not exist: {config['fantasia']['hf_cache_dir']}\n"
            "FANTASIA-Lite runs in HF_HUB_OFFLINE mode, so the ProtT5 weights "
            "must already be on disk before the rule is invoked."
        )
    if not config['fantasia']['lookup_dir']:
        raise ValueError(
            "fantasia.enable=1 but fantasia.lookup_dir is empty. Download the "
            "FANTASIA V1 lookup bundle from Zenodo record 17720428 and set "
            "fantasia.lookup_dir to the extracted directory in config.ini, or "
            "via GALBA2_FANTASIA_LOOKUP_DIR."
        )
    if not os.path.isdir(config['fantasia']['lookup_dir']):
        raise FileNotFoundError(
            f"fantasia.enable=1 but the configured lookup_dir does not exist: "
            f"{config['fantasia']['lookup_dir']}\n"
            "Download the FANTASIA V1 lookup bundle from "
            "https://zenodo.org/records/17720428/files/fantasia_lite_data_folder.zip "
            "and extract it to that path."
        )

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
# compleasm uses a flat layout ({path}/{lineage}/) while BUSCO uses
# {busco_download_path}/lineages/{lineage}/. By defaulting
# compleasm_download_path to the BUSCO "lineages" subdirectory, a single
# extracted lineage tarball serves both tools and we avoid downloading the
# same data twice. Users can still override this in [paths] if they want
# the caches kept apart.
config['compleasm_download_path'] = os.path.abspath(
    config_parser.get('paths', 'compleasm_download_path',
                      fallback=os.path.join(config['busco_download_path'], 'lineages'))
)

# OMAmer database path (for optional OMArk QC)
def _find_omamer_db():
    """Find LUCA.h5: [OMARK] omamer_db, shared_data, sibling repos."""
    explicit = config_parser.get('OMARK', 'omamer_db', fallback='') or \
               config_parser.get('PARAMS', 'omamer_db', fallback='')
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
config['rfam_path'] = os.path.abspath(
    config_parser.get('paths', 'rfam_path',
                      fallback=os.path.join(_shared_data, 'rfam'))
)

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
    if config['masking_tool'] == 'red':
        print("  ✓ Red (fast repeat detection and masking)")
    else:
        print("  ✓ RepeatModeler2 + RepeatMasker (genome masking)")
if config.get('use_minisplice', False):
    print("  ✓ Minisplice splice site scoring (CNN)")
print("  ✓ Miniprot protein-to-genome alignment")
print("  ✓ Miniprothint training gene and hint generation")
print("  ✓ AUGUSTUS training on protein-derived genes")
print("  ✓ AUGUSTUS prediction with protein hints")
if not config.get('disable_diamond_filter', False):
    print("  ✓ DIAMOND filtering of predictions against input proteins")
if config.get('run_ncrna', False):
    print("  ✓ ncRNA annotation (barrnap, tRNAscan-SE, Infernal/Rfam)")
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
    if config['masking_tool'] == 'red':
        include: "rules/preprocessing/run_red_masking.smk"
    else:
        include: "rules/preprocessing/run_masking.smk"

# Minisplice splice site scoring (optional)
if config.get('use_minisplice', False):
    include: "rules/miniprot/run_minisplice.smk"

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

# ncRNA annotation (optional)
RUN_NCRNA = config.get('run_ncrna', False)
if RUN_NCRNA:
    include: "rules/ncrna/run_barrnap.smk"
    include: "rules/ncrna/run_trnascan.smk"
    include: "rules/ncrna/run_infernal.smk"

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
