#!/usr/bin/env python3

import sys
import os
import re
import argparse
import shutil
import subprocess
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_transcript_filter import filter_transcripts  # noqa: E402

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023-2026. All rights reserved."
__credits__ = ""
__license__ = "Artistic License 2.0"
__version__ = "1.2.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

argparser = argparse.ArgumentParser(description = 'Filter genes predicted by ' +
                                    'AUGUSTUS within e.g. GALBA by simply ' +
                                    'discarding all genes that are not ' +
                                    'supported by a DIAMOND hit to the ' +
                                    'reference protein set. Warning: this is ' +
                                    'a poor approach if you do not have a ' +
                                    'high quality reference protein set from ' +
                                    'close relatives. This approach may lead ' +
                                    'to discarding species-specific gene ' +
                                    'models.')
argparser.add_argument('-r', '--ref_prots', type=str, required = True,
                            help = 'FASTA file with reference protien sequences')
argparser.add_argument('-g', '--gtf', type=str, required = False,
                            help = 'gtf file with gene models to be filtered')
argparser.add_argument('-a', '--pred_prots', type=str, required = False,
                            help = 'FASTA file with predicted protein sequences')
argparser.add_argument('-o', '--out_gtf', type=str, required = True,
                            help = 'Output file in gtf format')
argparser.add_argument('-t', '--threads', type=str, required = False, default = 1,
                       help = 'Number of threads to use for running DIAMOND')
argparser.add_argument('-d', '--diamond_dir', type=str, required = True,
                            help = 'Directory that contains DIAMOND binaries')

TX_RE = re.compile(r'transcript_id "([^"]+)"')


def clean_fasta(in_fasta, out_fasta):
    """
    Copy a protein FASTA file without blank lines and carriage returns and
    return the number of sequences.

    DIAMOND refuses a FASTA file whose first line is not a header ("Error:
    First line must begin with '>'"), e.g. OrthoDB partitions that start
    with an empty line. Exits if the file has no sequences or text before
    the first header.
    """
    n_seqs = 0
    try:
        with open(in_fasta, "r") as fin, open(out_fasta, "w") as fout:
            for line in fin:
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                if line.startswith(">"):
                    n_seqs += 1
                elif n_seqs == 0:
                    print("ERROR: " + in_fasta + " has sequence data before the "
                          + "first FASTA header: " + line[:60], file=sys.stderr)
                    sys.exit(1)
                fout.write(line + "\n")
    except IOError as e:
        print("ERROR: protein file could not be read: " + str(e), file=sys.stderr)
        sys.exit(1)
    if n_seqs == 0:
        print("ERROR: " + in_fasta + " contains no FASTA sequences", file=sys.stderr)
        sys.exit(1)
    return n_seqs


def run_diamond(ref_prots, pred_prots, threads, diamond_dir, work_dir):
    """
    Run DIAMOND to search for homologous proteins in the reference protein
    set. All files go to work_dir, so parallel jobs started from the same
    directory do not overwrite each other's results. Exits if DIAMOND fails;
    an empty hit table from a failed run would discard every gene.
    """
    ref_clean = os.path.join(work_dir, "ref_prots.clean.fa")
    n_ref = clean_fasta(ref_prots, ref_clean)
    print("Reference proteins: " + str(n_ref) + " sequences in " + ref_prots)
    diamond_out_file = os.path.join(work_dir, "diamond_out.tsv")
    diamond_cmd = [os.path.join(diamond_dir, "diamond"), "blastp",
                   "-d", ref_clean, "-q", pred_prots, "-o", diamond_out_file,
                   "-p", str(threads),
                   "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                   "evalue", "bitscore"]
    print("Running DIAMOND with command: " + " ".join(diamond_cmd))
    sys.stdout.flush()
    try:
        result = subprocess.run(diamond_cmd)
    except OSError as e:
        print("ERROR: DIAMOND could not be started: " + str(e), file=sys.stderr)
        sys.exit(1)
    if result.returncode != 0:
        print("ERROR: DIAMOND failed with exit code " + str(result.returncode)
              + "; see the DIAMOND messages above", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(diamond_out_file):
        print("ERROR: DIAMOND finished but wrote no output file "
              + diamond_out_file, file=sys.stderr)
        sys.exit(1)
    return diamond_out_file

def filter_gtf(gtf, diamond_out_file, out_gtf):
    """
    Read DIAMOND output, identify transcript names that have at least one hit,
    read gtf file, discard all transcripts that do not have a hit and all
    genes left without transcripts, write filtered gtf file.

    Gene lines are re-spanned to the transcripts that remain (BRAKER4 #62:
    version 1.0.0 kept the gene line of a gene whose widest transcript was
    discarded, so galba.gtf had genes wider than their transcripts). Gene
    IDs come from the GTF, so IDs renamed by normalize_gtf.py (g1_2) work.
    """
    # reading DIAMOND output
    try:
        with open(diamond_out_file, "r") as diamond_out:
            diamond_hits_tx = set()
            for line in diamond_out:
                fields = line.rstrip().split("\t")
                if len(fields) < 2 or fields[0].startswith("Running DIAMOND"):
                    continue
                diamond_hits_tx.add(fields[0])
    except IOError:
        print("ERROR: DIAMOND output file " + diamond_out_file + " could not be opened")
        sys.exit(1)

    # reading gtf file to filter out transcripts/genes that do not have a hit
    try:
        dropped_tx, dropped_genes = filter_transcripts(gtf, out_gtf, keep=diamond_hits_tx)
    except IOError as e:
        print("ERROR: gtf file could not be processed: " + str(e))
        sys.exit(1)
    for txid in dropped_tx:
        print("Discarding " + txid)
    print("Discarded " + str(len(dropped_tx)) + " transcripts without DIAMOND hit and "
          + str(len(dropped_genes)) + " genes left without transcripts")


def count_transcripts(gtf):
    """Number of distinct transcript IDs in a GTF file."""
    txs = set()
    with open(gtf, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or line.startswith("#"):
                continue
            m = TX_RE.search(fields[8])
            if m:
                txs.add(m.group(1))
            elif fields[2] == "transcript":
                txs.add(fields[8].strip())
    return len(txs)


def main():
    args = argparser.parse_args()
    out_dir = os.path.dirname(os.path.abspath(args.out_gtf))
    work_dir = tempfile.mkdtemp(prefix="diamond_filter_", dir=out_dir)
    try:
        diamond_out_file = run_diamond(args.ref_prots, args.pred_prots,
                                       args.threads, args.diamond_dir, work_dir)
        filter_gtf(args.gtf, diamond_out_file, args.out_gtf)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
    n_in = count_transcripts(args.gtf)
    n_out = count_transcripts(args.out_gtf)
    print("Transcripts before DIAMOND filter: " + str(n_in) + ", after: " + str(n_out))
    if n_in > 0 and n_out == 0:
        print("ERROR: none of the " + str(n_in) + " predicted transcripts has a "
              + "DIAMOND hit to " + args.ref_prots + ". Check that this is the "
              + "right protein file, or set disable_diamond_filter to keep the "
              + "AUGUSTUS predictions unfiltered.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
