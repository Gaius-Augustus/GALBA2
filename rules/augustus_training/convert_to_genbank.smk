rule convert_to_genbank:
    """
    Convert miniprot training genes GTF to GenBank format for AUGUSTUS training.

    In GALBA2, training genes come directly from miniprothint.py (the miniprot
    protein alignment pipeline). Unlike BRAKER, there is no GeneMark filtering
    step — the training genes are protein-derived and already high quality.

    Uses gff2gbSmallDNA.pl to build GenBank entries with flanking regions
    around each gene.

    Input:
        gtf: Miniprot training genes GTF
        flanking_value: Flanking region size (in bp) from compute_flanking_region
        genome: Genome assembly FASTA file

    Output:
        gb: GenBank file containing training genes
    """
    input:
        gtf = "output/{sample}/miniprot/miniprot_trainingGenes.gtf",
        flanking_value = "output/{sample}/flanking_dna_value.txt",
        genome = lambda w: os.path.join(get_output_dir(w), "genome.fa")
    output:
        gb = "output/{sample}/bonafide.gb"
    benchmark:
        "benchmarks/{sample}/convert_to_genbank/convert_to_genbank.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GALBA_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # Read the flanking DNA value from the file
        FLANKING=$(cat {input.flanking_value})
        echo "[INFO] Using flanking DNA value: $FLANKING"

        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.gb})

        # Convert training genes GTF to GenBank format
        gff2gbSmallDNA.pl {input.gtf} {input.genome} $FLANKING {output.gb}
        N_GB=$(grep -c '^LOCUS' {output.gb} || echo 0)
        echo "[INFO] gff2gbSmallDNA.pl produced $N_GB GenBank entries"

        if [ $N_GB -eq 0 ]; then
            echo "[ERROR] No training genes converted to GenBank format."
            echo "        This may be caused by:"
            echo "        (a) Too distant protein sequences for miniprot alignment"
            echo "        (b) Complex FASTA headers in the genome file"
            exit 1
        fi

        echo "[INFO] GTF to GenBank conversion completed"
        echo "[INFO] Output written to {output.gb}"
        """
