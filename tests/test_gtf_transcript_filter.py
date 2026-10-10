"""
Tests for scripts/gtf_transcript_filter.py and the two CLIs using it
(filter_gtf_by_txid.py, filter_gtf_by_diamond_against_ref.py), BRAKER4 #62.

After dropping transcripts, transcript lines must span their features, gene
lines must span their remaining transcripts, and genes without transcripts
must disappear. All test data is synthetic and generated in-memory.
"""

import importlib.util
import os
import subprocess
import sys

import pytest

SCRIPTS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
CHECK = os.path.join(SCRIPTS, "check_gtf_loci.py")
TXID_CLI = os.path.join(SCRIPTS, "filter_gtf_by_txid.py")
sys.path.insert(0, SCRIPTS)

from gtf_transcript_filter import filter_transcripts  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _line(seq, feature, start, end, strand, attr, source="AUGUSTUS"):
    frame = "0" if feature == "CDS" else "."
    return "\t".join([seq, source, feature, str(start), str(end), ".",
                      strand, frame, attr]) + "\n"


def _feat(seq, feature, start, end, strand, tx, gene):
    return _line(seq, feature, start, end, strand,
                 f'transcript_id "{tx}"; gene_id "{gene}";')


def _bare_tx(gene, tx, start, end, seq="X1", strand="+"):
    """Transcript with two exons, bare gene-less transcript line."""
    mid = start + (end - start) // 3
    return (_line(seq, "transcript", start, end, strand, tx)
            + _feat(seq, "exon", start, mid, strand, tx, gene)
            + _feat(seq, "CDS", start, mid, strand, tx, gene)
            + _feat(seq, "exon", end - 100, end, strand, tx, gene)
            + _feat(seq, "CDS", end - 100, end, strand, tx, gene))


def _attr_tx(gene, tx, start, end, seq="X1", strand="+"):
    return (_line(seq, "transcript", start, end, strand,
                  f'transcript_id "{tx}"; gene_id "{gene}";')
            + _feat(seq, "exon", start, end, strand, tx, gene)
            + _feat(seq, "CDS", start, end, strand, tx, gene))


def _bare_gene(gene, start, end, seq="X1", strand="+"):
    return _line(seq, "gene", start, end, strand, gene)


def _attr_gene(gene, start, end, seq="X1", strand="+"):
    return _line(seq, "gene", start, end, strand, f'gene_id "{gene}";')


HEADER = "# comment header\n# second line\n"

# g1: two transcripts, t1 wide (1000-5000), t2 narrow (2000-3000)
G1_BARE = (_bare_gene("g1", 1000, 5000) + _bare_tx("g1", "g1.t1", 1000, 5000)
           + _bare_tx("g1", "g1.t2", 2000, 3000))
# g2: single transcript, must stay untouched
G2_BARE = _bare_gene("g2", 8000, 9000) + _bare_tx("g2", "g2.t1", 8000, 9000)

G1_ATTR = (_attr_gene("g1", 1000, 5000) + _attr_tx("g1", "g1.t1", 1000, 5000)
           + _attr_tx("g1", "g1.t2", 2000, 3000))
G2_ATTR = _attr_gene("g2", 8000, 9000) + _attr_tx("g2", "g2.t1", 8000, 9000)


def _rows(path, feature=None):
    out = []
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if feature is None or f[2] == feature:
            out.append(f)
    return out


def _span_of(path, feature, ident):
    for f in _rows(path, feature):
        if f[8] == ident or f'"{ident}"' in f[8]:
            return int(f[3]), int(f[4])
    raise AssertionError(f"{feature} {ident} not found")


def _run_filter(tmp_path, text, **kw):
    gtf = tmp_path / "in.gtf"
    gtf.write_text(text)
    out = tmp_path / "out.gtf"
    res = filter_transcripts(str(gtf), str(out), **kw)
    return out, res


def _check(path):
    return subprocess.run([sys.executable, CHECK, str(path)],
                          capture_output=True, text=True)


@pytest.mark.parametrize("g1,g2", [(G1_BARE, G2_BARE), (G1_ATTR, G2_ATTR)],
                         ids=["bare", "attribute"])
def test_dropping_widest_transcript_shrinks_gene(tmp_path, g1, g2):
    out, (dtx, dgenes) = _run_filter(tmp_path, g1 + g2, drop={"g1.t1"})
    assert dtx == ["g1.t1"] and dgenes == []
    assert _span_of(out, "gene", "g1") == (2000, 3000)
    assert _span_of(out, "transcript", "g1.t2") == (2000, 3000)
    assert _span_of(out, "gene", "g2") == (8000, 9000)
    assert "g1.t1" not in out.read_text()
    assert _check(out).returncode == 0


def test_dropping_all_transcripts_removes_gene_others_byte_identical(tmp_path):
    out, (dtx, dgenes) = _run_filter(tmp_path, G1_BARE + G2_BARE,
                                     drop={"g1.t1", "g1.t2"})
    assert dtx == ["g1.t1", "g1.t2"] and dgenes == ["g1"]
    assert out.read_text() == G2_BARE


def test_keep_and_drop_semantics(tmp_path):
    text = G1_BARE + G2_BARE
    out_k, _ = _run_filter(tmp_path, text, keep={"g2.t1"})
    assert out_k.read_text() == G2_BARE
    out_d, (dtx, _) = _run_filter(tmp_path, text, drop={"g1.t1", "g1.t2"})
    assert out_d.read_text() == G2_BARE
    # empty keep drops everything, empty drop drops nothing
    out_e, (_, dg) = _run_filter(tmp_path, text, keep=set())
    assert out_e.read_text() == "" and dg == ["g1", "g2"]
    out_n, (dtx, dg) = _run_filter(tmp_path, text, drop=set())
    assert out_n.read_text() == text and dtx == [] and dg == []


def test_keep_and_drop_both_or_neither_raises(tmp_path):
    gtf = tmp_path / "in.gtf"
    gtf.write_text(G2_BARE)
    out = str(tmp_path / "out.gtf")
    with pytest.raises(ValueError):
        filter_transcripts(str(gtf), out, keep={"a"}, drop={"b"})
    with pytest.raises(ValueError):
        filter_transcripts(str(gtf), out)


def test_renamed_ids_attributed_to_renamed_gene(tmp_path):
    # g1_2 as produced by normalize_gtf.py; dropping g1_2.t1 must not
    # touch gene g1 and must shrink gene g1_2.
    text = (G1_BARE
            + _bare_gene("g1_2", 20000, 24000)
            + _bare_tx("g1_2", "g1_2.t1", 20000, 24000)
            + _bare_tx("g1_2", "g1_2.t2", 21000, 22000))
    out, (dtx, dgenes) = _run_filter(tmp_path, text, drop={"g1_2.t1"})
    assert dtx == ["g1_2.t1"] and dgenes == []
    assert _span_of(out, "gene", "g1_2") == (21000, 22000)
    assert _span_of(out, "gene", "g1") == (1000, 5000)
    assert _check(out).returncode == 0
    # dropping the only transcript of g1_2 removes g1_2, not g1
    text2 = (G1_BARE + _bare_gene("g1_2", 20000, 24000)
             + _bare_tx("g1_2", "g1_2.t1", 20000, 24000))
    out2, (_, dgenes2) = _run_filter(tmp_path, text2, drop={"g1_2.t1"})
    assert dgenes2 == ["g1_2"]
    assert out2.read_text() == G1_BARE


def test_renamed_ids_without_gene_id_on_features(tmp_path):
    # features carry only transcript_id; transcript belongs to preceding gene line
    def tx(t, s, e):
        return (_line("X1", "transcript", s, e, "+", t)
                + _line("X1", "exon", s, e, "+", f'transcript_id "{t}";'))
    text = (_bare_gene("g1_2", 100, 900) + tx("g1_2.t1", 100, 900)
            + tx("g1_2.t2", 200, 300))
    out, (_, dgenes) = _run_filter(tmp_path, text, drop={"g1_2.t1"})
    assert dgenes == []
    assert _span_of(out, "gene", "g1_2") == (200, 300)


def test_header_and_feature_order_preserved(tmp_path):
    out, _ = _run_filter(tmp_path, HEADER + G1_BARE + G2_BARE, drop={"g1.t1"})
    lines = out.read_text().splitlines(keepends=True)
    assert "".join(lines[:2]) == HEADER
    # g1.t2 block appears verbatim (already tight) in original feature order
    assert "".join(lines[2:8]) == _bare_gene("g1", 2000, 3000) + _bare_tx("g1", "g1.t2", 2000, 3000)
    assert [f[2] for f in _rows(out) if "g1.t2" in f[8]] == \
        ["transcript", "exon", "CDS", "exon", "CDS"]
    assert lines[8:] == (G2_BARE).splitlines(keepends=True)


def test_transcript_line_respanned_to_features(tmp_path):
    # transcript line wider than its features gets tightened
    text = (_bare_gene("g1", 1000, 5000)
            + _line("X1", "transcript", 1000, 5000, "+", "g1.t1")
            + _feat("X1", "exon", 1200, 1500, "+", "g1.t1", "g1")
            + _feat("X1", "exon", 4000, 4800, "+", "g1.t1", "g1"))
    out, _ = _run_filter(tmp_path, text, drop=set())
    assert _span_of(out, "transcript", "g1.t1") == (1200, 4800)
    assert _span_of(out, "gene", "g1") == (1200, 4800)


def test_txid_cli_end_to_end_then_loci_check(tmp_path):
    gtf = tmp_path / "in.gtf"
    gtf.write_text(HEADER + G1_BARE + G2_BARE
                   + _bare_gene("g1_2", 20000, 24000)
                   + _bare_tx("g1_2", "g1_2.t1", 20000, 24000)
                   + _bare_tx("g1_2", "g1_2.t2", 21000, 22000)
                   + _bare_gene("g3", 30000, 31000)
                   + _bare_tx("g3", "g3.t1", 30000, 31000))
    lst = tmp_path / "drop.lst"
    lst.write_text("g1.t1\ng1_2.t1\ng3.t1\n\n")
    out = tmp_path / "out.gtf"
    res = subprocess.run([sys.executable, TXID_CLI, "-g", str(gtf), "-l", str(lst),
                          "-o", str(out)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert _span_of(out, "gene", "g1") == (2000, 3000)
    assert _span_of(out, "gene", "g1_2") == (21000, 22000)
    ids = {f[8] for f in _rows(out, "gene")}
    assert ids == {"g1", "g2", "g1_2"}
    assert out.read_text().startswith(HEADER)
    chk = _check(out)
    assert chk.returncode == 0, chk.stderr
    # the unfiltered input is consistent too; sanity that checker sees this format
    assert _check(gtf).returncode == 0


def test_diamond_filter_gtf_keeps_hit_transcripts(tmp_path):
    mod = _load("filter_gtf_by_diamond_against_ref")
    gtf = tmp_path / "in.gtf"
    gtf.write_text(HEADER + G1_ATTR + G2_ATTR
                   + _attr_gene("g1_2", 20000, 24000)
                   + _attr_tx("g1_2", "g1_2.t1", 20000, 24000)
                   + _attr_tx("g1_2", "g1_2.t2", 21000, 22000))
    hits = tmp_path / "diamond.tsv"
    row = "\t".join(["{q}", "ref1", "90.0", "100", "5", "0", "1", "100", "1",
                     "100", "1e-30", "200.0"]) + "\n"
    # multiple hits for one query, plus an unrelated query not in the GTF
    hits.write_text("Running DIAMOND with command: fake\n"
                    + row.format(q="g1.t2") + row.format(q="g1.t2")
                    + row.format(q="g1_2.t2") + row.format(q="unrelated.t1"))
    out = tmp_path / "out.gtf"
    mod.filter_gtf(str(gtf), str(hits), str(out))

    txs = {f[8] for f in _rows(out, "transcript")}
    assert txs == {'transcript_id "g1.t2"; gene_id "g1";',
                   'transcript_id "g1_2.t2"; gene_id "g1_2";'}
    assert _span_of(out, "gene", "g1") == (2000, 3000)
    assert _span_of(out, "gene", "g1_2") == (21000, 22000)
    assert [f[8] for f in _rows(out, "gene")] == ['gene_id "g1";', 'gene_id "g1_2";']  # g2 dropped
    assert out.read_text().startswith(HEADER)
    chk = _check(out)
    assert chk.returncode == 0, chk.stderr


def test_diamond_filter_gtf_no_hits_gives_empty_output(tmp_path):
    mod = _load("filter_gtf_by_diamond_against_ref")
    gtf = tmp_path / "in.gtf"
    gtf.write_text(G1_BARE + G2_BARE)
    hits = tmp_path / "diamond.tsv"
    hits.write_text("")
    out = tmp_path / "out.gtf"
    mod.filter_gtf(str(gtf), str(hits), str(out))
    assert out.read_text() == ""


def test_diamond_filter_gtf_missing_hits_file_exits(tmp_path):
    mod = _load("filter_gtf_by_diamond_against_ref")
    gtf = tmp_path / "in.gtf"
    gtf.write_text(G2_BARE)
    with pytest.raises(SystemExit):
        mod.filter_gtf(str(gtf), str(tmp_path / "nope.tsv"), str(tmp_path / "o.gtf"))
