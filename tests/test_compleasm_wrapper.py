"""Tests for scripts/compleasm_wrapper.py: the placement-file marker that keeps
compleasm 0.2.7 from downloading (and crashing on) the --autolineage placement
files when a lineage is given explicitly. All data is synthetic."""

import os
import sys

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS)

from compleasm_wrapper import ensure_placement_marker, filter_long_proteins  # noqa: E402


def test_marker_created_for_explicit_lineage(tmp_path):
    lib = tmp_path / "lineages"
    lib.mkdir()
    marker = ensure_placement_marker(["run", "-l", "brassicales", "--odb", "odb12",
                                      "-a", "g.fa", "-L", str(lib)])
    assert marker == str(lib / "placement_files.done")
    assert (lib / "placement_files.done").is_file()
    assert (lib / "placement_files").is_dir()
    # idempotent
    assert ensure_placement_marker(["protein", "--library_path", str(lib)]) == marker


def test_no_marker_with_autolineage_or_without_library(tmp_path):
    lib = tmp_path / "lineages"
    lib.mkdir()
    assert ensure_placement_marker(["run", "--autolineage", "-L", str(lib)]) is None
    assert not (lib / "placement_files.done").exists()
    assert ensure_placement_marker(["run", "-l", "x_odb12"]) is None


def test_filter_long_proteins_drops_only_overlong(tmp_path):
    fa = tmp_path / "in.faa"
    fa.write_text(">short desc\nMKV*\n>huge\n" + "A" * 100000 + "\n>ok\nMKVL\n")
    out = tmp_path / "out.faa"
    dropped = filter_long_proteins(str(fa), str(out))
    assert dropped == ["huge"]
    assert out.read_text() == ">short desc\nMKV*\n>ok\nMKVL\n"
