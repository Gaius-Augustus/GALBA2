rule sanity_check_augustus:
    """
    Run sanity checks on AUGUSTUS output and create the galba.gtf file.

    In GALBA2, there is no TSEBRA merge step. The AUGUSTUS hints prediction
    output goes through:
    1. GTF sanity check (remove transcripts with overlapping CDS or CDS on
       different strands)
    2. Optional DIAMOND filtering against input proteins
    3. Rename to galba.gtf

    Input:
        augustus_gtf: AUGUSTUS predictions in GTF format
        proteins: Protein sequences for DIAMOND filtering (if enabled)

    Output:
        galba_gtf: Final gene predictions (galba.gtf)
    """
    input:
        augustus_gtf = "output/{sample}/galba.normalized.gtf",
        genome = lambda w: os.path.join(get_output_dir(w), "genome.fa"),
        proteins = lambda w: get_protein_fasta(w.sample)
    output:
        galba_gtf = "output/{sample}/galba.gtf"
    benchmark:
        "benchmarks/{sample}/sanity_check_augustus/sanity_check_augustus.txt"
    params:
        output_dir = lambda w: get_output_dir(w),
        disable_diamond_filter = config.get('disable_diamond_filter', False),
        translation_table = config.get("translation_table", 1)
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== POSTPROCESSING AUGUSTUS OUTPUT ====="

        # Step 1: Sanity check - remove transcripts with overlapping CDS features
        # or CDS features on different strands within one transcript
        echo "[INFO] Running GTF sanity check..."
        python3 $(which gtf_sanity_check.py) \
            -f {input.augustus_gtf} \
            -o {params.output_dir}/bad_transcripts.lst \
            > {params.output_dir}/gtf_sanity_check.log \
            2> {params.output_dir}/gtf_sanity_check.err

        BAD_COUNT=0
        if [ -f {params.output_dir}/bad_transcripts.lst ]; then
            BAD_COUNT=$(wc -l < {params.output_dir}/bad_transcripts.lst || echo 0)
        fi
        echo "[INFO] Found $BAD_COUNT problematic transcripts"

        if [ $BAD_COUNT -gt 0 ]; then
            python3 $(which filter_gtf_by_txid.py) \
                -g {input.augustus_gtf} \
                -l {params.output_dir}/bad_transcripts.lst \
                -o {params.output_dir}/augustus.clean.gtf \
                > {params.output_dir}/filter_gtf_by_txid.log \
                2> {params.output_dir}/filter_gtf_by_txid.err
        else
            cp {input.augustus_gtf} {params.output_dir}/augustus.clean.gtf
        fi

        # Step 2: Optional DIAMOND filter
        if [ "{params.disable_diamond_filter}" = "False" ]; then
            echo "[INFO] Running DIAMOND filter against input proteins..."

            # Extract protein sequences from AUGUSTUS predictions
            STEM="{params.output_dir}/augustus.hints"
            python3 {script_dir}/getAnnoFastaFromJoingenes.py \
                -g {input.genome} \
                -f {params.output_dir}/augustus.clean.gtf \
                -o $STEM \
                -t {params.translation_table} \
                -d {params.output_dir} \
                1> {params.output_dir}/getAnnoFasta_hints.stdout \
                2> {params.output_dir}/getAnnoFasta_hints.stderr

            # Run DIAMOND filtering
            python3 $(which filter_gtf_by_diamond_against_ref.py) \
                -r {input.proteins} \
                -g {params.output_dir}/augustus.clean.gtf \
                -o {output.galba_gtf} \
                -a {params.output_dir}/augustus.hints.aa \
                -t {threads} \
                -d /opt/diamond \
                1> {params.output_dir}/diamond_filter.stdout \
                2> {params.output_dir}/diamond_filter.stderr

            # Cleanup temporary files
            rm -f {params.output_dir}/augustus.hints.aa \
                  {params.output_dir}/augustus.hints.codingseq

            BEFORE=$(grep -cP '\tgene\t' {params.output_dir}/augustus.clean.gtf || echo 0)
            AFTER=$(grep -cP '\tgene\t' {output.galba_gtf} || echo 0)
            echo "[INFO] DIAMOND filter: $BEFORE genes -> $AFTER genes"

            # Report
            REPORT_DIR=output/{wildcards.sample}
            source {script_dir}/report_citations.sh
            cite diamond "$REPORT_DIR"
        else
            echo "[INFO] DIAMOND filter disabled, using AUGUSTUS output directly"
            mv {params.output_dir}/augustus.clean.gtf {output.galba_gtf}
        fi

        rm -f {params.output_dir}/augustus.clean.gtf

        GENES=$(grep -cP '\tgene\t' {output.galba_gtf} || echo 0)
        echo "[INFO] Final gene count: $GENES"
        echo "[INFO] Output: {output.galba_gtf}"
        echo "[INFO] ======================================="
        """


rule extract_final_sequences:
    """
    Extract coding sequences and proteins from final galba.gtf predictions.

    Input:
        galba_gtf: Final gene predictions GTF
        genome: Genome FASTA file

    Output:
        galba_codingseq: Final coding sequences
        galba_aa: Final protein sequences
    """
    input:
        galba_gtf = "output/{sample}/galba.gtf",
        genome = lambda w: os.path.join(get_output_dir(w), "genome.fa")
    output:
        galba_codingseq = "output/{sample}/galba.codingseq",
        galba_aa = "output/{sample}/galba.aa"
    benchmark:
        "benchmarks/{sample}/extract_final_sequences/extract_final_sequences.txt"
    params:
        output_dir = lambda w: get_output_dir(w),
        translation_table = config.get("translation_table", 1)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== EXTRACTING FINAL SEQUENCES ====="

        STEM=$(echo {output.galba_codingseq} | sed 's/\.codingseq$//')
        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {input.galba_gtf} \
            -o $STEM \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1> {params.output_dir}/getAnnoFasta_final.stdout \
            2> {params.output_dir}/getAnnoFasta_final.stderr

        GENES_GTF=$(grep -cP '\tgene\t' {input.galba_gtf} || echo 0)
        CODING_SEQ=$(grep -c "^>" {output.galba_codingseq} || echo 0)
        AA_SEQ=$(grep -c "^>" {output.galba_aa} || echo 0)

        echo "[INFO] Final gene predictions:"
        echo "[INFO]   Genes: $GENES_GTF"
        echo "[INFO]   Coding sequences: $CODING_SEQ"
        echo "[INFO]   Protein sequences: $AA_SEQ"

        # Verify no internal in-frame stop codons
        awk '
        /^>/ {{
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{ internal_stops++ }}
            seq=""
        }}
        !/^>/ {{seq=seq$0}}
        END {{
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{ internal_stops++ }}
            print internal_stops + 0
        }}
        ' {output.galba_aa} > {params.output_dir}/internal_stops_verify.txt

        INTERNAL_STOPS=$(cat {params.output_dir}/internal_stops_verify.txt)
        if [ "$INTERNAL_STOPS" -gt 0 ]; then
            echo "[ERROR] Found $INTERNAL_STOPS sequences with internal stop codons!"
            exit 1
        else
            echo "[INFO] Verification passed: No internal stop codons found"
            rm -f {params.output_dir}/internal_stops_verify.txt
        fi

        echo "[INFO] ======================================="
        """


rule assess_completeness:
    """
    Assess gene set completeness using compleasm (BUSCO assessment).
    """
    input:
        galba_aa = "output/{sample}/galba.aa"
    output:
        compleasm_summary = "output/{sample}/compleasm_proteins/summary.txt",
        compleasm_log = "output/{sample}/compleasm_proteins.log"
    benchmark:
        "benchmarks/{sample}/assess_completeness/assess_completeness.txt"
    params:
        busco_lineage = lambda w: get_busco_lineage(w),
        compleasm_outdir = lambda w: f"output/{w.sample}/compleasm_proteins",
        library_path = config['compleasm_download_path']
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== ASSESSING GENE SET COMPLETENESS WITH COMPLEASM =====" | tee {output.compleasm_log}
        echo "[INFO] BUSCO lineage: {params.busco_lineage}" | tee -a {output.compleasm_log}

        PROTEIN_COUNT=$(grep -c "^>" {input.galba_aa} || echo 0)
        echo "[INFO] Total proteins: $PROTEIN_COUNT" | tee -a {output.compleasm_log}

        if [ $PROTEIN_COUNT -eq 0 ]; then
            echo "[WARNING] No proteins found - skipping compleasm" | tee -a {output.compleasm_log}
            mkdir -p {params.compleasm_outdir}
            echo "No proteins available for assessment" > {output.compleasm_summary}
            exit 0
        fi

        if [ -d "{params.compleasm_outdir}" ]; then
            rm -rf {params.compleasm_outdir}
        fi

        mkdir -p {params.compleasm_outdir}
        mkdir -p {params.library_path}

        COMPLEASM_LINEAGE=$(echo "{params.busco_lineage}" | sed 's/_odb[0-9]*$/_odb12/')
        echo "[INFO] Running compleasm in protein mode..." | tee -a {output.compleasm_log}

        compleasm.py protein \
            -p {input.galba_aa} \
            -l $COMPLEASM_LINEAGE \
            -t {threads} \
            -o {params.compleasm_outdir} \
            -L {params.library_path} \
            2>&1 | tee -a {output.compleasm_log} || true

        if [ ! -f {output.compleasm_summary} ]; then
            FOUND_SUMMARY=$(find {params.compleasm_outdir} -name "summary.txt" 2>/dev/null | head -1)
            if [ -n "$FOUND_SUMMARY" ] && [ -f "$FOUND_SUMMARY" ]; then
                cp "$FOUND_SUMMARY" {output.compleasm_summary}
            else
                mkdir -p {params.compleasm_outdir}
                echo "Compleasm assessment failed" > {output.compleasm_summary}
            fi
        fi

        echo "[INFO] =======================================" | tee -a {output.compleasm_log}
        """
