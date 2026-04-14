"""
Optional minisplice splice site scoring for improved miniprot alignment.

minisplice (Yang, Huang & Li, 2025) uses a 1D CNN to score canonical
splice sites (GT donors, AG acceptors) on the genome. The scores are
passed to miniprot via --spsc to improve protein-to-genome alignment,
especially for distant homologs.

Only included when use_minisplice = 1 in config.ini.
Requires miniprot >= 0.14 and a pre-trained minisplice model.

Container: katharinahoff/galba2-tools:latest (must include minisplice)
"""


rule run_minisplice:
    """Score canonical splice sites on the genome using minisplice CNN."""
    input:
        genome = lambda w: get_masked_genome(w.sample)
    output:
        scores = "output/{sample}/minisplice/splice_scores.tsv"
    benchmark:
        "benchmarks/{sample}/minisplice/minisplice.txt"
    params:
        model = config.get('minisplice_model', '')
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_TOOLS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.scores})

        echo "[INFO] ===== MINISPLICE: SCORING SPLICE SITES ====="

        # Find minisplice binary
        MINISPLICE=$(which minisplice 2>/dev/null || true)
        if [ -z "$MINISPLICE" ]; then
            echo "[ERROR] minisplice not found in PATH"
            echo "[ERROR] Install minisplice (https://github.com/lh3/minisplice)"
            exit 1
        fi

        # Find model files (vi2-7k.kan + vi2-7k.kan.cali)
        MODEL_DIR="{params.model}"
        if [ -z "$MODEL_DIR" ] || [ ! -f "$MODEL_DIR" ]; then
            # Search default locations for vi2-7k.kan
            for candidate_dir in \
                {script_dir}/../shared_data/minisplice \
                /usr/local/share/minisplice; do
                if [ -f "$candidate_dir/vi2-7k.kan" ]; then
                    MODEL_DIR="$candidate_dir"
                    break
                fi
            done
        fi

        MODEL_KAN="$MODEL_DIR/vi2-7k.kan"
        MODEL_CALI="$MODEL_DIR/vi2-7k.kan.cali"

        if [ ! -f "$MODEL_KAN" ]; then
            echo "[ERROR] minisplice model vi2-7k.kan not found"
            echo "[ERROR] Download with: wget -O- https://zenodo.org/records/15931054/files/vi2-7k.tgz | tar zxf -"
            exit 1
        fi

        echo "[INFO] Model: $MODEL_KAN"
        echo "[INFO] Calibration: $MODEL_CALI"
        echo "[INFO] Threads: {threads}"

        # Run minisplice predict (-c for calibration file)
        CALI_ARG=""
        if [ -f "$MODEL_CALI" ]; then
            CALI_ARG="-c $MODEL_CALI"
        fi
        $MINISPLICE predict -t{threads} $CALI_ARG "$MODEL_KAN" {input.genome} > {output.scores}

        N_SCORES=$(wc -l < {output.scores})
        echo "[INFO] Scored $N_SCORES splice sites"

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        MINISPLICE_VER=$($MINISPLICE version 2>&1 | head -1 || echo "unknown")
        mkdir -p $(dirname "$VERSIONS_FILE")
        ( flock 9; printf "minisplice\t%s\n" "$MINISPLICE_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite minisplice "$REPORT_DIR"
        """
