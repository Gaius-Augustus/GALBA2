#!/usr/bin/env python3
"""
Generate GALBA hintsfile from miniprothint output.

Combines three sources:
1. CDSpart hints from training gene CDS features (padded by 15bp, src=C)
2. High-confidence (HC) hints from hc.gff (src=C, chain-grouped)
3. Lower-confidence (LC) hints from miniprothint.gff (src=P, with multiplicity)
"""

import sys
import re
import argparse


def main():
    parser = argparse.ArgumentParser(description="Generate GALBA hintsfile")
    parser.add_argument("--training_genes", required=True, help="Training genes GTF")
    parser.add_argument("--hc_gff", required=True, help="High-confidence hints GFF")
    parser.add_argument("--lc_gff", required=True, help="Lower-confidence hints GFF")
    parser.add_argument("--output", required=True, help="Output hintsfile")
    args = parser.parse_args()

    hints = []

    # Step 1: CDSpart hints from training genes (padded by 15bp)
    cds_count = 0
    with open(args.training_genes) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            # Extract transcript_id
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not m:
                continue
            grp = m.group(1)
            start = int(fields[3]) + 15
            end = int(fields[4]) - 15
            if start > end:
                start = end = (start + end) // 2
            hints.append(f"{fields[0]}\t{fields[1]}\tCDSpart\t{start}\t{end}\t{fields[5]}\t{fields[6]}\t{fields[7]}\tgrp={grp};src=C;pri=4;")
            cds_count += 1
    print(f"[INFO] CDSpart hints from training genes: {cds_count}", file=sys.stderr)

    # Step 2: High-confidence hints from hc.gff
    grps = {}
    with open(args.hc_gff) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Replace start_codon/stop_codon
            line = line.replace('\tstart_codon\t', '\tstart\t')
            line = line.replace('\tstop_codon\t', '\tstop\t')
            m = re.match(r'(.*)\tal_score=\d+\.\d+;.*(prots=[^;]*)', line)
            if m:
                first_part = m.group(1)
                prot_key = m.group(2)
                grps.setdefault(prot_key, []).append(first_part)

    grp_id = 1
    hc_count = 0
    for key, lines in grps.items():
        for l in lines:
            hints.append(f"{l}\tgrp=chain_{grp_id};src=C;pri=4;")
            hc_count += 1
        grp_id += 1
    print(f"[INFO] High-confidence hints: {hc_count}", file=sys.stderr)

    # Step 3: Lower-confidence hints from miniprothint.gff
    lc_count = 0
    with open(args.lc_gff) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = line.replace('\tstart_codon\t', '\tstart\t')
            line = line.replace('\tstop_codon\t', '\tstop\t')
            m = re.match(r'(.*)\t.*prots=([^;]+);', line)
            if m:
                first_part = m.group(1)
                prots = m.group(2)
                mult = prots.count(',') + 1
                hints.append(f"{first_part}\tmult={mult};src=P;pri=4;")
                lc_count += 1
    print(f"[INFO] Lower-confidence hints: {lc_count}", file=sys.stderr)

    # Sort hints
    def sort_key(h):
        fields = h.split('\t')
        return (fields[0], fields[2], int(fields[3]), int(fields[4]))

    hints.sort(key=sort_key)

    # Write output
    with open(args.output, 'w') as f:
        for h in hints:
            f.write(h + '\n')

    print(f"[INFO] Total hints written: {len(hints)}", file=sys.stderr)


if __name__ == '__main__':
    main()
