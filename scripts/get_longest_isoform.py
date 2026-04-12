#!/usr/bin/env python3
"""
Select the longest coding isoform per gene locus from a GTF file.

For each gene, keeps only the transcript with the longest total CDS length.
This avoids inflated BUSCO duplicate counts from alternative splicing.

Usage:
    python3 get_longest_isoform.py -g input.gtf -o output.gtf
"""

import argparse
import sys
import re
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Select longest coding isoform per gene locus")
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output")
    args = parser.parse_args()

    # Parse GTF: collect CDS lengths per transcript, and gene-to-transcript mapping
    gene_to_tx = defaultdict(list)  # gene_id -> [tx_id, ...]
    tx_cds_len = defaultdict(int)   # tx_id -> total CDS length
    tx_to_gene = {}                 # tx_id -> gene_id

    lines_by_tx = defaultdict(list)  # tx_id -> [lines]
    gene_lines = defaultdict(list)   # gene_id -> [gene-level lines]

    with open(args.gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = fields[8]

            gene_match = re.search(r'gene_id "([^"]+)"', attrs)
            tx_match = re.search(r'transcript_id "([^"]+)"', attrs)

            if not gene_match:
                continue

            gene_id = gene_match.group(1)

            if feature == 'gene':
                gene_lines[gene_id].append(line)
                continue

            if not tx_match:
                continue

            tx_id = tx_match.group(1)
            tx_to_gene[tx_id] = gene_id

            if tx_id not in gene_to_tx[gene_id]:
                gene_to_tx[gene_id].append(tx_id)

            lines_by_tx[tx_id].append(line)

            if feature == 'CDS':
                start = int(fields[3])
                end = int(fields[4])
                tx_cds_len[tx_id] += end - start + 1

    # For each gene, select the transcript with the longest CDS
    kept_tx = set()
    for gene_id, tx_list in gene_to_tx.items():
        if not tx_list:
            continue
        best_tx = max(tx_list, key=lambda t: tx_cds_len.get(t, 0))
        kept_tx.add(best_tx)

    # Write output
    with open(args.output, 'w') as out:
        for gene_id in gene_to_tx:
            # Write gene line
            for line in gene_lines.get(gene_id, []):
                out.write(line)
            # Write lines for the kept transcript only
            for tx_id in gene_to_tx[gene_id]:
                if tx_id in kept_tx:
                    for line in lines_by_tx[tx_id]:
                        out.write(line)

    if not args.quiet:
        n_genes = len(gene_to_tx)
        n_tx_in = sum(len(v) for v in gene_to_tx.values())
        n_tx_out = len(kept_tx)
        print(f"Genes: {n_genes}, Transcripts in: {n_tx_in}, Transcripts out: {n_tx_out}",
              file=sys.stderr)


if __name__ == '__main__':
    main()
