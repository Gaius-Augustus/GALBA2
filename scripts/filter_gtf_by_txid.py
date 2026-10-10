#!/usr/bin/env python3

__author__ = "Katharina J. Hoff with a little help of ChatGPT"
__copyright__ = "Copyright 2022-2026. All rights reserved."
__credits__ = "Katharina Hoff"
__version__ = "1.1.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

"""
Remove the transcripts listed in a text file from an AUGUSTUS GTF.

Gene lines are re-spanned to the transcripts that remain, and genes left
without transcripts are dropped (see gtf_transcript_filter.py). Version
1.0.1 copied gene lines through unchanged, so a gene whose widest
transcript was removed kept a span larger than its transcripts in
galba.gtf.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_transcript_filter import filter_transcripts  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description='Remove transcripts from a gtf file')
    parser.add_argument('-g', '--gtf', required=True, type=str, help="Name of (sorted) gtf file")
    parser.add_argument('-l', '--exclude_lst', required=True, type=str, help="Text file with transcript IDs")
    parser.add_argument('-o', '--output_file', required=True, type=str, help="Output file in gtf format")
    args = parser.parse_args()

    with open(args.exclude_lst) as f:
        transcripts_to_remove = {line.strip() for line in f if line.strip()}

    dropped_tx, dropped_genes = filter_transcripts(
        args.gtf, args.output_file, drop=transcripts_to_remove)
    print(f"Removed {len(dropped_tx)} of {len(transcripts_to_remove)} listed "
          f"transcripts and {len(dropped_genes)} genes left without transcripts")


if __name__ == "__main__":
    main()
