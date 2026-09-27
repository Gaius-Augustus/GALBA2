"""
barrnap rRNA gene prediction and merging with GALBA annotations.

Runs barrnap on the genome to predict ribosomal RNA genes (18S, 28S, 5.8S, 5S),
then merges the rRNA features (along with tRNA, Infernal, and lncRNA results
when also enabled) into the final GFF3 annotation.

Only included when run_ncrna = 1 in config.ini.

Container: quay.io/biocontainers/barrnap:0.9--hdfd78af_4
"""


rule run_barrnap:
    """Run barrnap to predict rRNA genes on the genome."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        gff="output/{sample}/barrnap/rRNA.gff3"
    log:
        "logs/{sample}/barrnap/barrnap.log"
    benchmark:
        "benchmarks/{sample}/barrnap/barrnap.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BARRNAP_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff})

        barrnap --kingdom euk --threads {threads} \
            {input.genome} > {output.gff}.tmp 2> {log} || true

        # Prefix rRNA gene IDs with sample name and add unique counter
        if [ -s {output.gff}.tmp ] && grep -qv '^#' {output.gff}.tmp; then
            # POSIX-portable Name extraction: the biocontainer may ship busybox
            # awk, which does NOT support gawk's 3-arg match($str, regex, array).
            # Use 2-arg match() with RSTART/RLENGTH + substr() instead.
            awk -F'\t' -v OFS='\t' -v p="{wildcards.sample}" '
                BEGIN {{n=1}}
                /^#/ {{print; next}}
                {{
                    oldname = ""
                    if (match($9, /Name=[^;]+/)) {{
                        oldname = substr($9, RSTART + 5, RLENGTH - 5)
                    }}
                    newname = p "-rRNA_" n "_" oldname
                    gsub(/Name=[^;]+/, "Name=" newname, $9)
                    if ($9 ~ /ID=/) {{
                        gsub(/ID=[^;]+/, "ID=" newname, $9)
                    }} else {{
                        $9 = "ID=" newname ";" $9
                    }}
                    n++
                    print
                }}
            ' {output.gff}.tmp > {output.gff}
        else
            echo "##gff-version 3" > {output.gff}
            echo "No rRNA genes found" >> {log}
        fi
        rm -f {output.gff}.tmp

        n_rrna=$(grep -cv '^#' {output.gff} || echo 0)
        echo "barrnap predicted $n_rrna rRNA features" >> {log}

        # Record software version (LC_ALL=C avoids locale warnings in biocontainer)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        # BusyBox grep in this container has no -P; use sed for portable extraction
        BARRNAP_VER=$(LC_ALL=C barrnap --version 2>&1 | sed -n 's/.*\([0-9][0-9]*\.[0-9][0-9]*[^ ]*\).*/\1/p' | head -1 || true)
        ( flock 9; printf "barrnap\t%s\n" "$BARRNAP_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite barrnap "$REPORT_DIR" || true
        """


def _get_merge_ncrna_inputs(wildcards):
    """Get inputs for merging rRNA, tRNA, Infernal, and lncRNA into GFF3."""
    sample = wildcards.sample
    inputs = {
        "gff3": f"output/{sample}/galba.gff3",
        "rrna": f"output/{sample}/barrnap/rRNA.gff3",
        "trna": f"output/{sample}/ncrna/tRNAs.gff3",
        "infernal": f"output/{sample}/ncrna/ncRNAs_infernal.gff3",
    }
    # lncRNA (FEELnc) requires transcript evidence — not available in GALBA2.
    return inputs


rule merge_rrna_into_gff3:
    """Merge barrnap rRNA, tRNA, Infernal, and lncRNA into the final GFF3."""
    input:
        unpack(_get_merge_ncrna_inputs)
    output:
        merged="output/{sample}/galba_with_ncRNA.gff3"
    log:
        "logs/{sample}/barrnap/merge_rrna.log"
    benchmark:
        "benchmarks/{sample}/merge_rrna_into_gff3/merge_rrna_into_gff3.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        # highest priority first; lncRNA only when the sample has one
        ncrna=lambda wildcards, input: [input.rrna, input.trna, input.infernal]
            + ([input.lncrna] if "lncrna" in input.keys() else [])
    shell:
        """
        set -euo pipefail
        # Protein-coding genes unchanged, then ncRNA genes in priority order
        # (rRNA > tRNA > Infernal > lncRNA). Each ncRNA becomes gene -> RNA ->
        # exon with gene_biotype on the gene (NCBI/Ensembl GFF3, Annotrieve).
        # ncRNA genes overlapping coding exons (>50%, same strand) or a
        # higher-priority ncRNA gene (>50%, same strand) are dropped.
        python3 {script_dir}/merge_ncrna_gff3.py \
            --coding {input.gff3} \
            --ncrna {params.ncrna} \
            -o {output.merged} \
            2> {log}
        """
