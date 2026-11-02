from genefamily.hmmer import parse_tblout, select_hits


TBL = """# target name        accession  query name  accession  E-value  score  bias
G1.1                  -          Cyclin_N    PF00134.29 1e-40    150.0  0.1  desc
G1.2                  -          Cyclin_N    PF00134.29 1e-30    120.0  0.1  desc
G2                    -          Cyclin_N    PF00134.29 3e-12     60.0  0.2  desc
G3                    -          Cyclin_N    PF00134.29 0.4        5.0   0.3  desc
#\n
"""


def test_parse_tblout_skips_comments(tmp_path):
    path = tmp_path / "h.tblout"
    path.write_text(TBL)
    hits = parse_tblout(str(path))
    assert len(hits) == 4
    assert hits[0]["target"] == "G1.1"
    assert hits[0]["evalue"] == 1e-40


def test_select_hits_applies_evalue_and_isoform_collapse(tmp_path):
    path = tmp_path / "h.tblout"
    path.write_text(TBL)
    hits = parse_tblout(str(path))
    selected = select_hits(hits, max_evalue=1e-5, best_per_gene=True)
    names = [h["target"] for h in selected]
    assert names == ["G1.1", "G2"]          # G1.2 collapsed, G3 filtered


def test_select_hits_can_keep_isoforms(tmp_path):
    path = tmp_path / "h.tblout"
    path.write_text(TBL)
    hits = parse_tblout(str(path))
    selected = select_hits(hits, max_evalue=1e-5, best_per_gene=False)
    assert len(selected) == 3
