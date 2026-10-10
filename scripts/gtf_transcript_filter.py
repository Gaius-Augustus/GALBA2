#!/usr/bin/env python3

"""
Drop transcripts from an AUGUSTUS-style GTF and keep gene lines consistent.

Shared by filter_gtf_by_txid.py and filter_gtf_by_diamond_against_ref.py.
Both used to copy gene lines through unchanged when they removed some of
a gene's transcripts, so the gene line kept the span of the dropped
transcript (BRAKER4 #62 symptom in galba.gtf, found in an A. thaliana run:
3 of 30549 genes wider than their transcripts). check_gtf_loci.py only ran
on galba.normalized.gtf, upstream of both filters, and could not see it.

After filtering, every transcript line spans exactly its features and every
gene line spans exactly its remaining transcripts. Genes without transcripts
are dropped. Gene and transcript IDs are taken from the GTF itself, never
derived from the ID pattern, so gene IDs renamed by normalize_gtf.py
(<gene_id>_2, ...) work too.

Handles attribute-style column 9 as well as the bare AUGUSTUS gene and
transcript lines, where column 9 holds only the ID. Output keeps the input
order of genes, transcripts and features, grouped per gene.
"""

import re
from collections import OrderedDict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')
TX_SUFFIX_RE = re.compile(r'\.t\d+$')


def _attr(attr, pattern):
    m = pattern.search(attr)
    return m.group(1) if m else None


def parse_gtf(gtf_file):
    """Return (genes, header) where genes maps gene_id to
    {"line": fields or None, "transcripts": OrderedDict(tx_id -> {"line":
    fields or None, "features": [fields, ...]})}. header holds comment
    lines verbatim."""
    header = []
    gene_lines = {}
    tx_lines = {}
    features = OrderedDict()
    tx_gene = {}
    tx_order = []
    gene_order = []
    last_gene = None

    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature == 'gene':
                gene_id = _attr(fields[8], GENE_RE) or fields[8].strip()
                gene_lines[gene_id] = fields
                if gene_id not in gene_order:
                    gene_order.append(gene_id)
                last_gene = gene_id
            elif feature in ('transcript', 'mRNA'):
                tx_id = _attr(fields[8], TX_RE) or fields[8].strip()
                tx_lines[tx_id] = fields
                if tx_id not in features:
                    features[tx_id] = []
                    tx_order.append(tx_id)
                gene_id = _attr(fields[8], GENE_RE)
                if gene_id:
                    tx_gene[tx_id] = gene_id
                elif last_gene is not None and tx_id not in tx_gene:
                    # Bare AUGUSTUS transcript line: belongs to the gene
                    # line written just before it. A later feature line
                    # with gene_id overrides this.
                    tx_gene[tx_id] = last_gene
            else:
                tx_id = _attr(fields[8], TX_RE)
                if tx_id is None:
                    continue
                if tx_id not in features:
                    features[tx_id] = []
                    tx_order.append(tx_id)
                features[tx_id].append(fields)
                gene_id = _attr(fields[8], GENE_RE)
                if gene_id:
                    tx_gene[tx_id] = gene_id

    genes = OrderedDict()
    for gene_id in gene_order:
        genes[gene_id] = {"line": gene_lines[gene_id],
                          "transcripts": OrderedDict()}
    for tx_id in tx_order:
        gene_id = tx_gene.get(tx_id) or TX_SUFFIX_RE.sub('', tx_id)
        if gene_id not in genes:
            genes[gene_id] = {"line": None, "transcripts": OrderedDict()}
        genes[gene_id]["transcripts"][tx_id] = {
            "line": tx_lines.get(tx_id), "features": features[tx_id]}
    return genes, header


def _span(rows):
    return (min(int(r[3]) for r in rows), max(int(r[4]) for r in rows))


def _with_span(fields, start, end):
    out = list(fields)
    out[3], out[4] = str(start), str(end)
    return out


def filter_transcripts(gtf_in, gtf_out, keep=None, drop=None):
    """Write gtf_in to gtf_out without the transcripts in `drop`, or without
    the transcripts not in `keep` (exactly one of the two must be given).

    Transcript and gene lines are re-spanned to their remaining children;
    genes left without transcripts are removed.

    Returns (dropped_transcripts, dropped_genes) as lists of IDs in GTF order.
    """
    if (keep is None) == (drop is None):
        raise ValueError("give exactly one of keep or drop")
    genes, header = parse_gtf(gtf_in)
    dropped_tx, dropped_genes = [], []

    with open(gtf_out, "w") as out:
        out.writelines(header)
        for gene_id, gene in genes.items():
            kept = OrderedDict()
            for tx_id, tx in gene["transcripts"].items():
                remove = (tx_id in drop) if drop is not None else (tx_id not in keep)
                if remove:
                    dropped_tx.append(tx_id)
                else:
                    kept[tx_id] = tx
            if not kept:
                dropped_genes.append(gene_id)
                continue

            tx_rows = {}
            for tx_id, tx in kept.items():
                if tx["features"]:
                    start, end = _span(tx["features"])
                    line = tx["line"]
                    tx_rows[tx_id] = (_with_span(line, start, end)
                                      if line is not None else None, start, end)
                elif tx["line"] is not None:
                    tx_rows[tx_id] = (tx["line"], int(tx["line"][3]), int(tx["line"][4]))
                else:
                    tx_rows[tx_id] = (None, None, None)

            if gene["line"] is not None:
                spans = [(s, e) for _, s, e in tx_rows.values() if s is not None]
                if spans:
                    out.write("\t".join(_with_span(gene["line"],
                                                   min(s for s, _ in spans),
                                                   max(e for _, e in spans))) + "\n")
                else:
                    out.write("\t".join(gene["line"]) + "\n")
            for tx_id, tx in kept.items():
                line = tx_rows[tx_id][0]
                if line is not None:
                    out.write("\t".join(line) + "\n")
                for row in tx["features"]:
                    out.write("\t".join(row) + "\n")

    return dropped_tx, dropped_genes
