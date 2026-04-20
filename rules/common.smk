"""
Common functions and configuration for GALBA2 Snakemake workflow.

This module provides:
- Sample parsing from samples.csv
- Input validation (protein evidence required)
- Helper functions for dynamic input resolution

GALBA2 is a protein-only gene prediction pipeline. Every sample must
provide a genome and a protein FASTA file.
"""

import pandas as pd
import os
from pathlib import Path

# ==============================================================================
# Configuration and Samples
# ==============================================================================

# Read samples.csv
samples_csv = config.get("samples_file", "samples.csv")
if not os.path.isfile(samples_csv):
    raise FileNotFoundError(
        f"Samples file '{samples_csv}' not found. "
        "Please create it or set samples_file in config.ini [paths]. "
        "See the README for the expected CSV format."
    )
samples_df = pd.read_csv(samples_csv)

# Get list of all samples
SAMPLES = samples_df["sample_name"].tolist()

# ==============================================================================
# Data Type Detection
# ==============================================================================

def detect_data_types(sample):
    """
    Detect which data types are present for a sample.

    In GALBA2, proteins are required; the only optional flags are
    whether a pre-masked genome is provided and whether a reference GTF
    exists for evaluation.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]

    return {
        "has_proteins": pd.notna(row.get("protein_fasta")),
        "has_masked_genome": pd.notna(row.get("genome_masked")),
        "needs_masking": pd.isna(row.get("genome_masked")),
        "has_reference_gtf": pd.notna(row.get("reference_gtf"))
    }

# ==============================================================================
# Input Validation
# ==============================================================================

def validate_samples():
    """Validate sample configurations before workflow execution."""
    for idx, row in samples_df.iterrows():
        sample = row["sample_name"]

        # Rule 1: Genome is required
        if pd.isna(row.get("genome")):
            raise ValueError(f"{sample}: genome column is required")

        # Rule 2: Proteins are required for GALBA2
        if pd.isna(row.get("protein_fasta")):
            raise ValueError(
                f"{sample}: protein_fasta column is required for GALBA2. "
                "GALBA2 is a protein-only gene prediction pipeline."
            )

        # Rule 3: busco_lineage is required
        if pd.isna(row.get("busco_lineage")) or not str(row.get("busco_lineage")).strip():
            raise ValueError(f"{sample}: busco_lineage column is required in samples.csv")

# Run validation at workflow start
validate_samples()

# ==============================================================================
# Helper Functions
# ==============================================================================

def get_masked_genome(sample):
    """Get the masked genome file path for a sample.

    Returns:
    - Pre-masked genome (cleaned): output/{sample}/genome_masked.fa
    - RepeatMasker output: output/{sample}/preprocessing/genome.fa.masked
    - Cleaned unmasked genome: output/{sample}/genome.fa
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("genome_masked")):
        return f"output/{sample}/genome_masked.fa"
    if GLOBAL_DATA_TYPES.get("needs_masking", False):
        return f"output/{sample}/preprocessing/genome.fa.masked"
    return f"output/{sample}/genome.fa"

def get_genome(sample):
    """Get the cleaned genome file path (headers simplified to accession only)."""
    return f"output/{sample}/genome.fa"

def get_raw_genome(sample):
    """Get the original user-provided genome path (before header cleaning)."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row["genome"]

def get_protein_fasta_files(sample):
    """Get list of protein FASTA file paths for a sample.

    Supports colon-separated paths for multiple files.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("protein_fasta")):
        return str(row["protein_fasta"]).split(":")
    return []

def get_protein_fasta(sample):
    """Get the protein FASTA path for pipeline consumption.

    If multiple files are specified (colon-separated), returns the path
    to the merged output. If single file, returns it directly.
    """
    files = get_protein_fasta_files(sample)
    if len(files) > 1:
        return f"output/{sample}/preprocessing/proteins_merged.fa"
    elif len(files) == 1:
        return files[0]
    return None

def get_reference_gtf(sample):
    """Get reference annotation GTF path for a sample, or None."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row.get("reference_gtf") if pd.notna(row.get("reference_gtf")) else None

def get_busco_lineage(sample):
    """Get BUSCO lineage for a sample — must be set in samples.csv."""
    sample = sample.sample if hasattr(sample, 'sample') else sample
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    lineage = row.get("busco_lineage")
    if not lineage or (isinstance(lineage, float) and pd.isna(lineage)):
        raise ValueError(
            f"Sample '{sample}' is missing required 'busco_lineage' column in samples.csv"
        )
    return lineage

def get_output_dir(wildcards):
    """Get the output directory for a sample."""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    return f"output/{sample}"

def get_species_name(wildcards):
    """Get the species name for AUGUSTUS (derived from sample name)."""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    # Replace special characters with underscores for AUGUSTUS compatibility
    return sample.replace("-", "_").replace(".", "_")

# ==============================================================================
# Container Configuration
# ==============================================================================

# Main containers:
# - GALBA_TOOLS: miniprot, boundary_scorer, miniprothint, compleasm (small, ~50 MB)
# - AUGUSTUS_CONTAINER: augustus, etraining, diamond, all Perl/Python scripts (249 MB)
# Until the galba2-tools container is built and pushed, the galba-notebook
# container can be used as a fallback for all rules.
GALBA_TOOLS_CONTAINER = config.get("galba_tools_image", "docker://katharinahoff/galba-notebook:latest")
AUGUSTUS_CONTAINER = config.get("augustus_image", "docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9716f88_9")
# Alias for rules that haven't been migrated yet
GALBA_CONTAINER = config.get("galba_image", "docker://katharinahoff/galba-notebook:latest")
GFFCOMPARE_CONTAINER = config.get("gffcompare_image", "docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1")
RED_CONTAINER = config.get("red_image", "docker://quay.io/biocontainers/red:2018.09.10--h9948957_3")
BARRNAP_CONTAINER = config.get("barrnap_image", "docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4")
AGAT_CONTAINER = config.get("agat_image", "docker://quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0")
BUSCO_CONTAINER = config.get("busco_image", "docker://ezlabgva/busco:v6.0.0_cv1")
OMARK_CONTAINER = config.get("omark_image", "docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0")
TETOOLS_CONTAINER = config.get("tetools_image", "docker://dfam/tetools:latest")

# ==============================================================================
# Workflow Summary
# ==============================================================================

def print_workflow_summary():
    """Print summary of detected samples and their configurations."""
    print("\n" + "="*70)
    print("GALBA2 Snakemake Workflow")
    print("="*70)
    print(f"\nTotal samples: {len(SAMPLES)}")

    for sample in SAMPLES:
        types = detect_data_types(sample)
        print(f"\n  {sample}:")
        print(f"    Mode: protein homology")
        print(f"    Genome: {'pre-masked' if types['has_masked_genome'] else 'unmasked'}")
        n_prots = len(get_protein_fasta_files(sample))
        print(f"    Proteins: {n_prots} file(s)")
        if types['has_reference_gtf']:
            print(f"    Reference GTF: yes (evaluation enabled)")

    print("\n" + "="*70 + "\n")

# Print summary when workflow loads
print_workflow_summary()
