"""
Merge and decompress protein FASTA files.

When multiple protein files are specified (colon-separated in samples.csv),
they must be concatenated before being passed to miniprot and the DIAMOND
filter. A single gzipped file (.gz, e.g. an OrthoDB partition as
downloaded) also goes through this rule, because the DIAMOND filter reads
only plain FASTA (ported from BRAKER4 d661e6b, Gaius-Augustus/BRAKER#789).

This is a concatenation without blank lines — no deduplication.
"""


rule merge_proteins:
    """Concatenate protein FASTA files into one, decompressing .gz files."""
    input:
        proteins=lambda wildcards: get_protein_fasta_files(wildcards.sample)
    output:
        merged="output/{sample}/preprocessing/proteins_merged.fa"
    log:
        "logs/{sample}/merge_proteins/merge.log"
    benchmark:
        "benchmarks/{sample}/merge_proteins/merge_proteins.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})
        # awk instead of cat: drops blank lines (DIAMOND rejects a FASTA
        # whose first line is not a header) and ends every file with a
        # newline, so a file without a final newline is not glued to the
        # next file's first header
        # gzip -cdf decompresses .gz files and passes plain files through
        # unchanged (gzipped protein_fasta, Gaius-Augustus/BRAKER#789)
        for f in {input.proteins}; do
            gzip -cdf -- "$f" | awk '{{ sub(/\r$/, "") }} NF'
        done > {output.merged}
        n_seqs=$(grep -c '^>' {output.merged})
        echo "Merged/decompressed $(echo {input.proteins} | wc -w) protein file(s): $n_seqs sequences" > {log}
        """
