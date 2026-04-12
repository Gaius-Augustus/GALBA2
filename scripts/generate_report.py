#!/usr/bin/env python3
"""
Generate a self-contained HTML report for a GALBA2 pipeline run.

Usage:
    python3 generate_report.py -d <output_dir> -o <results_dir> -s <sample_name>

Output:
  - galba_report.html: Self-contained HTML with inline CSS
  - galba_citations.bib: Deduplicated BibTeX file
"""

import argparse
import os
import re
import sys
from pathlib import Path


def read_file(path):
    """Read file content, return empty string if not found."""
    try:
        with open(path) as f:
            return f.read()
    except (FileNotFoundError, PermissionError):
        return ""


def count_features(gtf_path):
    """Count features in a GTF file."""
    counts = {}
    try:
        with open(gtf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    feat = fields[2]
                    counts[feat] = counts.get(feat, 0) + 1
    except FileNotFoundError:
        pass
    return counts


def count_seqs(fasta_path):
    """Count sequences in a FASTA file."""
    try:
        with open(fasta_path) as f:
            return sum(1 for line in f if line.startswith('>'))
    except FileNotFoundError:
        return 0


def dedup_bib(bib_path, out_path):
    """Deduplicate BibTeX entries by key."""
    entries = {}
    current_key = None
    current_lines = []
    try:
        with open(bib_path) as f:
            for line in f:
                m = re.match(r'@\w+\{(\w+),', line)
                if m:
                    if current_key and current_key not in entries:
                        entries[current_key] = ''.join(current_lines)
                    current_key = m.group(1)
                    current_lines = [line]
                elif current_key is not None:
                    current_lines.append(line)
            if current_key and current_key not in entries:
                entries[current_key] = ''.join(current_lines)
    except FileNotFoundError:
        pass
    with open(out_path, 'w') as f:
        for entry in entries.values():
            f.write(entry + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--outdir", required=True)
    parser.add_argument("-o", "--results", required=True)
    parser.add_argument("-s", "--sample", required=True)
    args = parser.parse_args()

    outdir = args.outdir
    results = args.results
    sample = args.sample

    # Gather stats
    gtf_path = os.path.join(outdir, "galba.gtf")
    features = count_features(gtf_path)
    n_genes = features.get("gene", 0)
    n_transcripts = features.get("transcript", 0)
    n_cds = features.get("CDS", 0)
    n_exons = features.get("exon", 0)
    n_proteins = count_seqs(os.path.join(outdir, "galba.aa"))
    n_coding = count_seqs(os.path.join(outdir, "galba.codingseq"))

    # Read compleasm summary
    compleasm_prot = read_file(os.path.join(outdir, "compleasm_proteins", "summary.txt"))
    compleasm_genome = read_file(os.path.join(outdir, "compleasm_genome_out", "summary.txt"))

    # Read citations
    citations_txt = read_file(os.path.join(outdir, "report_citations.txt"))

    # Read training stats
    etraining_count = read_file(os.path.join(outdir, "etraining_gene_count.txt")).strip()
    locus_count = read_file(os.path.join(outdir, "locus_count.txt")).strip()

    # Software versions
    versions = read_file(os.path.join(outdir, "software_versions.tsv"))

    # Build HTML
    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>GALBA2 Report &mdash; {sample}</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; }}
h1 {{ color: #2c5f2d; border-bottom: 2px solid #2c5f2d; padding-bottom: 10px; }}
h2 {{ color: #3a7d44; margin-top: 30px; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background: #f0f7f0; }}
.stat {{ font-size: 1.2em; font-weight: bold; color: #2c5f2d; }}
pre {{ background: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 0.9em; }}
.section {{ margin: 20px 0; padding: 15px; background: #fafafa; border-radius: 5px; }}
</style>
</head>
<body>
<h1>GALBA2 Report</h1>
<p>Sample: <strong>{sample}</strong></p>
<p>Pipeline: GALBA2 &mdash; protein homology based genome annotation</p>

<h2>Gene Prediction Summary</h2>
<table>
<tr><th>Metric</th><th>Count</th></tr>
<tr><td>Genes</td><td class="stat">{n_genes}</td></tr>
<tr><td>Transcripts</td><td>{n_transcripts}</td></tr>
<tr><td>CDS features</td><td>{n_cds}</td></tr>
<tr><td>Exon features</td><td>{n_exons}</td></tr>
<tr><td>Protein sequences</td><td>{n_proteins}</td></tr>
<tr><td>Coding sequences</td><td>{n_coding}</td></tr>
</table>

<h2>Training Statistics</h2>
<div class="section">
<p>Etraining gene filtering: {etraining_count if etraining_count else 'N/A'}</p>
<p>DIAMOND redundancy removal: {locus_count if locus_count else 'N/A'}</p>
</div>

<h2>Completeness Assessment</h2>
<h3>Genome (compleasm)</h3>
<pre>{compleasm_genome if compleasm_genome else 'Not available'}</pre>
<h3>Predicted Proteins (compleasm)</h3>
<pre>{compleasm_prot if compleasm_prot else 'Not available'}</pre>

<h2>Software Versions</h2>
<pre>{versions if versions else 'Not available'}</pre>

<h2>Citations</h2>
<pre>{citations_txt if citations_txt else 'Not available'}</pre>

<h2>Output Files</h2>
<ul>
<li><code>galba.gtf</code> &mdash; Gene predictions in GTF format</li>
<li><code>galba.gff3</code> &mdash; Gene predictions in GFF3 format</li>
<li><code>galba.aa</code> &mdash; Protein sequences</li>
<li><code>galba.codingseq</code> &mdash; Coding sequences</li>
<li><code>hintsfile.gff</code> &mdash; Protein-derived hints used for prediction</li>
</ul>

<p><em>Generated by GALBA2</em></p>
</body>
</html>"""

    # Write HTML
    html_path = os.path.join(results, "galba_report.html")
    with open(html_path, 'w') as f:
        f.write(html)
    print(f"[INFO] Report written to {html_path}")

    # Dedup BibTeX
    bib_in = os.path.join(outdir, "report_citations.bib")
    bib_out = os.path.join(results, "galba_citations.bib")
    dedup_bib(bib_in, bib_out)
    print(f"[INFO] BibTeX written to {bib_out}")


if __name__ == "__main__":
    main()
