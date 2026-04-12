rule generate_galba_hints:
    """
    Generate the final AUGUSTUS hintsfile from miniprothint output.

    This rule combines three sources of hints (matching original GALBA logic):
    1. CDSpart hints derived from training gene CDS features (padded by 15bp)
    2. High-confidence (HC) hints from miniprothint hc.gff
    3. Lower-confidence (LC) hints from miniprothint.gff

    Input:
        training_genes: Training genes GTF from miniprothint
        hc_gff: High-confidence hints from miniprothint
        lc_gff: Lower-confidence hints from miniprothint

    Output:
        hintsfile: Final hintsfile for AUGUSTUS prediction
    """
    input:
        training_genes = "output/{sample}/miniprot/miniprot_trainingGenes.gtf",
        hc_gff = "output/{sample}/miniprot/hc.gff",
        lc_gff = "output/{sample}/miniprot/miniprothint.gff"
    output:
        hintsfile = "output/{sample}/hintsfile.gff"
    benchmark:
        "benchmarks/{sample}/generate_galba_hints/generate_galba_hints.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== GENERATING GALBA HINTS ====="

        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        python3 {script_dir}/generate_galba_hints.py \
            --training_genes {input.training_genes} \
            --hc_gff {input.hc_gff} \
            --lc_gff {input.lc_gff} \
            --output {output.hintsfile}

        FINAL_HINTS=$(wc -l < {output.hintsfile})
        echo "[INFO] Final hintsfile: $FINAL_HINTS hints"
        echo "[INFO] Output: {output.hintsfile}"
        echo "[INFO] ======================================="
        """
