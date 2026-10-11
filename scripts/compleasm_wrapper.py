#!/usr/bin/env python3
"""
Run compleasm.py with two workarounds for compleasm 0.2.8.

1. Offline fallback (issue #105). compleasm defines its own
   `class URLError(OSError)`, which shadows urllib.error.URLError. Its
   `except URLError:` around the file_versions.tsv.hash download therefore
   never catches what urllib raises, and compleasm aborts on any network or
   SSL failure even when file_versions.tsv and the lineage are already in
   --library_path. The wrapper re-raises urllib's URLError as compleasm's own
   class, so compleasm prints a warning and uses the cached files.

2. Over-long proteins (issues #91, #92). hmmsearch aborts on any target
   sequence of 100,000 aa or more ("Target sequence length > 100K, over
   comparison pipeline limit"), which makes `compleasm protein` fail as a
   whole. In protein mode the wrapper passes compleasm a copy of the -p file
   without such sequences. No real protein comes close (titin is ~35,000 aa),
   so no BUSCO is lost.

3. Placement files (compleasm 0.2.7). Downloader.__init__ downloads the
   "placement files" of --autolineage unconditionally, and its parser crashes
   on the file_versions.tsv the BUSCO server has served since mid-2026
   ("ValueError: too many values to unpack" in download_placement), leaving a
   placement_files.tmp behind that makes every later run stop with "another
   process is downloading". compleasm skips the download when
   <library_path>/placement_files.done exists, so the wrapper creates that
   marker when a lineage is given (no --autolineage) and --library_path is
   set. An --autolineage run removes the marker itself and downloads again.

Usage:
    compleasm_wrapper.py [/path/to/compleasm.py] <compleasm arguments>

Without a path, compleasm.py is looked up in PATH and /opt/compleasm_kit.
"""

import os
import runpy
import shutil
import sys
import tempfile
import urllib.error
import urllib.request

# hmmsearch (p7_pipeline.c) refuses targets longer than 100,000 residues
MAX_PROTEIN_LEN = 100000


def filter_long_proteins(in_path, out_path, max_len=MAX_PROTEIN_LEN):
    """Copy FASTA in_path to out_path, dropping records with >= max_len residues.

    Returns the list of dropped record IDs.
    """
    dropped = []
    header = None
    seq = []

    def flush(out):
        if header is None:
            return
        s = "".join(seq)
        if len(s.rstrip("*")) >= max_len:
            dropped.append(header[1:].split()[0] if len(header) > 1 else "")
        else:
            out.write(header + "\n" + s + "\n")

    with open(in_path) as fin, open(out_path, "w") as out:
        for line in fin:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush(out)
                header = line
                seq = []
            elif line:
                seq.append(line.strip())
        flush(out)
    return dropped


def patch_urlretrieve():
    """Make urllib's URLError catchable by compleasm's `except URLError:`."""
    orig = urllib.request.urlretrieve

    def urlretrieve(*args, **kwargs):
        try:
            return orig(*args, **kwargs)
        except urllib.error.URLError as e:
            # runpy executes compleasm as the __main__ module, so its own
            # URLError class lives there.
            own = getattr(sys.modules.get("__main__"), "URLError", None)
            if own is None or own is urllib.error.URLError:
                raise
            raise own(str(e)) from e

    urllib.request.urlretrieve = urlretrieve


def ensure_placement_marker(cargs):
    """Create <library_path>/placement_files.done unless --autolineage is given.

    Returns the marker path when it was created or already present, else None.
    """
    if "--autolineage" in cargs:
        return None
    lib = None
    for flag in ("-L", "--library_path"):
        if flag in cargs:
            i = cargs.index(flag) + 1
            if i < len(cargs):
                lib = cargs[i]
            break
    if not lib:
        return None
    marker = os.path.join(lib, "placement_files.done")
    os.makedirs(os.path.join(lib, "placement_files"), exist_ok=True)
    if not os.path.exists(marker):
        open(marker, "w").close()
    return marker


def find_compleasm():
    for cand in (shutil.which("compleasm.py"), shutil.which("compleasm"),
                 "/opt/compleasm_kit/compleasm.py"):
        if cand and os.path.isfile(cand):
            return cand
    sys.exit("ERROR: compleasm.py not found in PATH or /opt/compleasm_kit")


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: compleasm_wrapper.py [/path/to/compleasm.py] <compleasm arguments>")
    if os.path.isfile(sys.argv[1]):
        compleasm = os.path.realpath(sys.argv[1])
        cargs = sys.argv[2:]
    else:
        compleasm = os.path.realpath(find_compleasm())
        cargs = sys.argv[1:]

    tmp_fasta = None
    if cargs and cargs[0] == "protein":
        for flag in ("-p", "--proteins"):
            if flag in cargs:
                i = cargs.index(flag) + 1
                if i < len(cargs) and os.path.isfile(cargs[i]):
                    fd, tmp_fasta = tempfile.mkstemp(prefix="compleasm_input_", suffix=".faa")
                    os.close(fd)
                    dropped = filter_long_proteins(cargs[i], tmp_fasta)
                    if dropped:
                        print(f"[compleasm_wrapper] Skipped {len(dropped)} protein(s) of "
                              f">= {MAX_PROTEIN_LEN} aa that hmmsearch cannot handle: "
                              + ", ".join(dropped), file=sys.stderr)
                        cargs[i] = tmp_fasta
                    else:
                        os.remove(tmp_fasta)
                        tmp_fasta = None
                break

    ensure_placement_marker(cargs)
    patch_urlretrieve()
    # compleasm imports _version from its own directory and locates miniprot
    # and hmmsearch relative to __file__, which run_path sets correctly.
    sys.path.insert(0, os.path.dirname(compleasm))
    sys.argv = [compleasm] + cargs
    try:
        runpy.run_path(compleasm, run_name="__main__")
    finally:
        if tmp_fasta and os.path.exists(tmp_fasta):
            os.remove(tmp_fasta)


if __name__ == "__main__":
    main()
