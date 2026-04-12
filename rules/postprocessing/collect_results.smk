"""
Collect final results into a clean directory and remove intermediate files.

The results directory contains everything a user needs:
- Gene predictions (GTF, GFF3, protein and CDS sequences)
- Quality control reports (BUSCO, compleasm, OMArk, gffcompare)
- Evidence support summary
"""


def _get_collect_inputs(wildcards):
    """Get all final output files that must exist before collection."""
    sample = wildcards.sample
    inputs = {
        "gtf": f"output/{sample}/galba.gtf",
        "aa": f"output/{sample}/galba.aa",
        "codingseq": f"output/{sample}/galba.codingseq",
        "gff3": f"output/{sample}/galba.gff3",
        # gene_support requires intervaltree which is not in the GALBA container
        # "gene_support": f"output/{sample}/gene_support.tsv",
        "headers_fixed": f"output/{sample}/preprocessing/.headers_fixed"
    }

    if not config.get('skip_busco', False):
        inputs["busco_summary"] = f"output/{sample}/busco/busco_summary.txt"
        inputs["busco_figure"] = f"output/{sample}/busco/busco_figure.png"

    inputs["compleasm"] = f"output/{sample}/compleasm_proteins/summary.txt"
    inputs["compleasm_genome"] = f"output/{sample}/compleasm_genome_out/summary.txt"

    if config.get('run_omark', False):
        inputs["omark"] = f"output/{sample}/omark/omark_summary.txt"

    types = detect_data_types(sample)
    if types.get('has_reference_gtf'):
        inputs["gffcompare"] = f"output/{sample}/gffcompare/gffcompare.stats"

    return inputs


rule collect_results:
    """Collect important output files and clean up intermediates."""
    input:
        unpack(_get_collect_inputs)
    output:
        done="output/{sample}/results/.done"
    benchmark:
        "benchmarks/{sample}/collect_results/collect_results.txt"
    params:
        sample="{sample}",
        outdir="output/{sample}",
        resultsdir="output/{sample}/results",
        no_cleanup=config.get('no_cleanup', False)
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    container:
        GALBA_CONTAINER
    shell:
        r"""
        set -euo pipefail

        RESULTS="{params.resultsdir}"
        OUTDIR="{params.outdir}"

        mkdir -p "$RESULTS"
        mkdir -p "$RESULTS/quality_control"

        echo "[INFO] Collecting results for sample {params.sample}..."

        # --- Copy genome if headers were cleaned ---
        if [ -f "$OUTDIR/preprocessing/.headers_fixed" ] && grep -q "yes" "$OUTDIR/preprocessing/.headers_fixed" 2>/dev/null; then
            if [ -f "$OUTDIR/genome.fa" ]; then
                cp "$OUTDIR/genome.fa" "$RESULTS/"
                echo "[INFO] Including cleaned genome in results"
            fi
        fi

        # --- Core gene predictions (copy first, gzip AFTER report generation) ---
        for f in galba.gtf galba.gff3 galba.aa galba.codingseq; do
            if [ -f "$OUTDIR/$f" ]; then
                cp "$OUTDIR/$f" "$RESULTS/"
            fi
        done

        # --- Evidence support ---
        cp "$OUTDIR/gene_support.tsv" "$RESULTS/" 2>/dev/null || true
        cp "$OUTDIR/hintsfile.gff"    "$RESULTS/" 2>/dev/null || true

        # --- Quality control ---

        # BUSCO
        if [ -f "$OUTDIR/busco/busco_summary.txt" ]; then
            cp "$OUTDIR/busco/busco_summary.txt" "$RESULTS/quality_control/"
        fi
        if [ -f "$OUTDIR/busco/busco_figure.png" ]; then
            cp "$OUTDIR/busco/busco_figure.png"  "$RESULTS/quality_control/"
        fi
        if [ -d "$OUTDIR/busco" ]; then
            for bmode in genome proteins; do
                if [ -d "$OUTDIR/busco/$bmode" ]; then
                    summary=$(find "$OUTDIR/busco/$bmode" -name "short_summary*.txt" 2>/dev/null | head -1 || true)
                    if [ -n "$summary" ] && [ -f "$summary" ]; then
                        cp "$summary" "$RESULTS/quality_control/busco_${{bmode}}_short_summary.txt"
                    fi
                fi
            done
        fi

        # Compleasm (proteome)
        if [ -f "$OUTDIR/compleasm_proteins/summary.txt" ]; then
            cp "$OUTDIR/compleasm_proteins/summary.txt" "$RESULTS/quality_control/compleasm_summary.txt"
        fi

        # Compleasm (genome)
        if [ -f "$OUTDIR/compleasm_genome_out/summary.txt" ]; then
            mkdir -p "$RESULTS/quality_control/compleasm_genome_out"
            cp "$OUTDIR/compleasm_genome_out/summary.txt" "$RESULTS/quality_control/compleasm_genome_out/summary.txt"
        fi

        # OMArk
        if [ -d "$OUTDIR/omark" ]; then
            if [ -f "$OUTDIR/omark/omark_summary.txt" ]; then
                cp "$OUTDIR/omark/omark_summary.txt" "$RESULTS/quality_control/"
            fi
            for f in "$OUTDIR"/omark/*_detailed_summary.txt "$OUTDIR"/omark/*.sum "$OUTDIR"/omark/*.tax; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/quality_control/"
                fi
            done
        fi

        # gffcompare
        if [ -d "$OUTDIR/gffcompare" ]; then
            cp "$OUTDIR/gffcompare/gffcompare.stats" "$RESULTS/quality_control/" 2>/dev/null || true
        fi

        # --- Software versions ---
        if [ -f "$OUTDIR/software_versions.tsv" ]; then
            awk -F'\t' '!seen[$1]++' "$OUTDIR/software_versions.tsv" | sort -f > "$RESULTS/software_versions.tsv"
        fi

        # --- Training summary ---
        echo "[INFO] Generating training summary..."
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        python3 {script_dir}/training_summary.py \
            -d "$OUTDIR" \
            -o "$RESULTS/quality_control" \
            2>/dev/null || echo "[WARNING] Training summary generation failed (non-fatal)"

        # --- Gene set statistics ---
        echo "[INFO] Generating gene set statistics..."
        SUPPORT_ARG=""
        if [ -f "$OUTDIR/gene_support.tsv" ]; then
            SUPPORT_ARG="-s $OUTDIR/gene_support.tsv"
        fi
        python3 {script_dir}/gene_set_statistics.py \
            -g "$OUTDIR/galba.gtf" \
            -o "$RESULTS/quality_control" \
            $SUPPORT_ARG \
            2>/dev/null || echo "[WARNING] Gene set statistics generation failed (non-fatal)"

        # --- Completeness plot ---
        python3 {script_dir}/completeness_plot.py \
            -d "$RESULTS/quality_control" \
            -o "$RESULTS/quality_control/completeness.png" \
            2>/dev/null || echo "[WARNING] Completeness plot generation failed (non-fatal)"

        # --- Generate HTML report ---
        python3 {script_dir}/generate_report.py \
            -d "$OUTDIR" \
            -o "$RESULTS" \
            -s "{params.sample}" \
            2>/dev/null || echo "[WARNING] HTML report generation failed (non-fatal)"

        # --- Gzip large result files ---
        echo "[INFO] Compressing result files..."
        for f in "$RESULTS"/galba.gtf "$RESULTS"/galba.gff3 \
                 "$RESULTS"/galba.aa "$RESULTS"/galba.codingseq \
                 "$RESULTS"/hintsfile.gff "$RESULTS"/genome.fa; do
            if [ -f "$f" ]; then
                gzip -f "$f"
            fi
        done

        # --- Clean up intermediate files ---
        if [ "{params.no_cleanup}" = "True" ]; then
            echo "[INFO] Skipping cleanup (no_cleanup=1 in config.ini)"
        else
            echo "[INFO] Removing intermediate files..."
            for item in "$OUTDIR"/*; do
                basename=$(basename "$item")
                if [ "$basename" = "results" ]; then
                    continue
                fi
                rm -rf "$item"
            done
            rm -f "$(pwd)"/*.agat.log "$(pwd)"/busco_*.log 2>/dev/null || true
        fi

        touch "$RESULTS/.done"

        echo "[INFO] Results collected in $RESULTS:"
        find "$RESULTS" -type f -not -name '.done' | sort | while read f; do
            size=$(du -h "$f" | cut -f1)
            rel=$(echo "$f" | sed "s|$RESULTS/||")
            echo "  $size  $rel"
        done
        """
