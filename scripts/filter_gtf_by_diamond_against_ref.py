#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_transcript_filter import filter_transcripts  # noqa: E402

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023-2026. All rights reserved."
__credits__ = ""
__license__ = "Artistic License 2.0"
__version__ = "1.1.0"
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


def run_diamond(ref_prots, pred_prots, threads, diamond_dir):
    """
    Run DIAMOND to search for homologous proteins in the reference protein set
    """
    diamond_out_file = "diamond_out.tsv"
    diamond_out = open(diamond_out_file, "w")
    diamond_cmd = diamond_dir + "/diamond blastp -d " + ref_prots + " -q " + pred_prots + " -o " + diamond_out_file + " -p " + threads + " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    diamond_out.write("Running DIAMOND with command: " + diamond_cmd + "\n")
    diamond_out.flush() 
    diamond_out.close()
    subprocess.call(diamond_cmd, shell=True)
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


def main():
    args = argparser.parse_args()
    diamond_out_file = run_diamond(args.ref_prots, args.pred_prots, args.threads, args.diamond_dir)
    filter_gtf(args.gtf, diamond_out_file, args.out_gtf)
    os.remove(diamond_out_file)

if __name__ == "__main__":
    main()
