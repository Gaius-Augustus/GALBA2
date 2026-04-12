rule miniprot_index:
    """
    Build a miniprot genome index for fast protein-to-genome alignment.

    Input:
        genome: Genome FASTA file (header-cleaned)

    Output:
        index: Miniprot genome index (.mpi)
    """
    input:
        genome = lambda w: os.path.join(get_output_dir(w), "genome.fa")
    output:
        index = "output/{sample}/miniprot/genome.mpi"
    benchmark:
        "benchmarks/{sample}/miniprot_index/miniprot_index.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_TOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== MINIPROT: BUILDING GENOME INDEX ====="
        echo "[INFO] Genome: {input.genome}"
        echo "[INFO] Threads: {threads}"

        mkdir -p $(dirname {output.index})

        miniprot -t{threads} -d {output.index} {input.genome} \
            > output/{wildcards.sample}/miniprot/miniprot_index.log \
            2> output/{wildcards.sample}/miniprot/miniprot_index.err

        echo "[INFO] Genome index created: {output.index}"

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        MINIPROT_VER=$(miniprot --version 2>&1 || true)
        mkdir -p $(dirname "$VERSIONS_FILE")
        ( flock 9; printf "miniprot\t%s\n" "$MINIPROT_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite miniprot "$REPORT_DIR"
        """


rule miniprot_align:
    """
    Align protein sequences to the genome using miniprot.

    Produces alignment output in miniprot's --aln format, which is then
    processed by miniprot_boundary_scorer and miniprothint.py.

    Input:
        index: Miniprot genome index
        proteins: Merged protein FASTA

    Output:
        aln: Miniprot alignment file (--aln format)
    """
    input:
        index = "output/{sample}/miniprot/genome.mpi",
        proteins = lambda w: get_protein_fasta(w.sample)
    output:
        aln = "output/{sample}/miniprot/miniprot.aln"
    benchmark:
        "benchmarks/{sample}/miniprot_align/miniprot_align.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_TOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== MINIPROT: ALIGNING PROTEINS ====="
        echo "[INFO] Index: {input.index}"
        echo "[INFO] Proteins: {input.proteins}"
        echo "[INFO] Threads: {threads}"

        miniprot -I -ut{threads} --outn=1 --aln \
            {input.index} {input.proteins} \
            > {output.aln} \
            2> output/{wildcards.sample}/miniprot/miniprot_align.err

        ALN_LINES=$(wc -l < {output.aln})
        echo "[INFO] Alignment produced $ALN_LINES lines"
        echo "[INFO] Alignment written to: {output.aln}"
        """


rule miniprot_boundary_scorer:
    """
    Score splice site boundaries using the miniprot_boundary_scorer tool.

    The boundary scorer uses a BLOSUM62 matrix to evaluate the quality of
    splice site boundaries in miniprot alignments, producing a scored GFF
    that is then processed by miniprothint.py.

    Input:
        aln: Miniprot alignment file

    Output:
        scored_gff: Scored GFF output from boundary scorer
    """
    input:
        aln = "output/{sample}/miniprot/miniprot.aln"
    output:
        scored_gff = "output/{sample}/miniprot/miniprot_scored.gff"
    benchmark:
        "benchmarks/{sample}/miniprot_boundary_scorer/miniprot_boundary_scorer.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_TOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== MINIPROT: BOUNDARY SCORING ====="

        # Find miniprot_boundary_scorer binary
        SCORER=$(which miniprot_boundary_scorer 2>/dev/null || true)
        if [ -z "$SCORER" ]; then
            echo "[ERROR] miniprot_boundary_scorer not found in PATH"
            exit 1
        fi

        # Find blosum62.csv — check BLOSUM62_PATH env var, then known locations
        BLOSUM="${{BLOSUM62_PATH:-}}"
        if [ -z "$BLOSUM" ] || [ ! -f "$BLOSUM" ]; then
            for candidate in \
                /usr/local/share/galba2/blosum62.csv \
                "$(dirname $SCORER)/blosum62.csv" \
                /opt/miniprot-boundary-scorer/blosum62.csv \
                /opt/GALBA/scripts/blosum62.csv; do
                if [ -f "$candidate" ]; then
                    BLOSUM="$candidate"
                    break
                fi
            done
        fi
        if [ ! -f "$BLOSUM" ]; then
            echo "[ERROR] blosum62.csv not found"
            exit 1
        fi

        echo "[INFO] Scorer: $SCORER"
        echo "[INFO] BLOSUM: $BLOSUM"

        $SCORER -o {output.scored_gff} -s $BLOSUM < {input.aln}

        GFF_LINES=$(wc -l < {output.scored_gff})
        echo "[INFO] Boundary scorer produced $GFF_LINES lines"
        """


rule run_miniprothint:
    """
    Run miniprothint.py to generate training genes and hints from miniprot output.

    miniprothint.py converts the boundary scorer output into:
    - miniprot_trainingGenes.gtf: Training gene structures for AUGUSTUS
    - hc.gff: High-confidence hints (chain-grouped)
    - miniprothint.gff: Lower-confidence hints with multiplicity

    Input:
        scored_gff: Scored GFF from miniprot_boundary_scorer

    Output:
        training_genes: Training genes GTF
        hc_gff: High-confidence hints
        lc_gff: Lower-confidence hints (miniprothint.gff)
    """
    input:
        scored_gff = "output/{sample}/miniprot/miniprot_scored.gff"
    output:
        training_genes = "output/{sample}/miniprot/miniprot_trainingGenes.gtf",
        hc_gff = "output/{sample}/miniprot/hc.gff",
        lc_gff = "output/{sample}/miniprot/miniprothint.gff"
    benchmark:
        "benchmarks/{sample}/run_miniprothint/run_miniprothint.txt"
    params:
        workdir = "output/{sample}/miniprot"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_TOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== MINIPROTHINT: GENERATING HINTS AND TRAINING GENES ====="

        # Find miniprothint.py
        MINIPROTHINT=""
        for dir in /opt/miniprot-boundary-scorer /opt/miniprothint /opt/ProtHint/dependencies /opt/GALBA/scripts /usr/local/bin; do
            if [ -f "$dir/miniprothint.py" ]; then
                MINIPROTHINT="$dir/miniprothint.py"
                break
            fi
        done

        if [ -z "$MINIPROTHINT" ]; then
            MINIPROTHINT=$(which miniprothint.py 2>/dev/null || true)
        fi

        if [ -z "$MINIPROTHINT" ]; then
            echo "[ERROR] miniprothint.py not found"
            exit 1
        fi

        echo "[INFO] miniprothint.py: $MINIPROTHINT"

        # --ignoreCoverage prints hints to hc.gff ignoring coverage if most hints have coverage = 1
        # --topNperSeed 10 and --minScoreFraction 0.5 match original GALBA defaults
        python3 $MINIPROTHINT {input.scored_gff} \
            --workdir {params.workdir} \
            --ignoreCoverage \
            --topNperSeed 10 \
            --minScoreFraction 0.5

        # Verify outputs
        TRAINING_GENES=$(grep -c $'\tCDS\t' {output.training_genes} || echo 0)
        HC_LINES=$(wc -l < {output.hc_gff})
        LC_LINES=$(wc -l < {output.lc_gff})

        echo "[INFO] Training genes CDS lines: $TRAINING_GENES"
        echo "[INFO] High-confidence hints: $HC_LINES lines"
        echo "[INFO] Lower-confidence hints: $LC_LINES lines"

        if [ "$TRAINING_GENES" -eq 0 ]; then
            echo "[ERROR] No training genes produced. Check protein alignment quality."
            exit 1
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite miniprothint "$REPORT_DIR"
        """
