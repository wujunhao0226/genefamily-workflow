from genefamily.domains import dedupe, parse_rpsblast, summarise


def _row(protein, acc, start, end, evalue="1e-20", title=None):
    title = title or f"{acc}, MODEL_{acc}, synthetic"
    return "\t".join([protein, f"gnl|CDD|{acc[2:]}", title, "30", "100", "5", "1",
                      str(start), str(end), "1", "100", evalue, "100"])


def test_parse_rpsblast_reads_identity_from_title(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(_row("P1", "cd00001", 5, 120) + "\n")
    hits = parse_rpsblast(str(path))
    assert len(hits) == 1
    assert hits[0]["accession"] == "cd00001"
    assert hits[0]["name"] == "MODEL_cd00001"
    assert hits[0]["superfamily"] is False


def test_dedupe_keeps_best_overlapping_specific_hit():
    hits = [
        {"protein": "P1", "accession": "cd1", "name": "A", "qstart": 10, "qend": 120,
         "evalue": 1e-30, "bitscore": 100, "superfamily": False},
        {"protein": "P1", "accession": "cd2", "name": "B", "qstart": 20, "qend": 110,
         "evalue": 1e-20, "bitscore": 80, "superfamily": False},
    ]
    kept = dedupe(hits)
    assert len(kept) == 1
    assert kept[0]["accession"] == "cd1"


def test_dedupe_drops_superfamily_overlapping_specific():
    hits = [
        {"protein": "P1", "accession": "cd1", "name": "A", "qstart": 10, "qend": 120,
         "evalue": 1e-30, "bitscore": 100, "superfamily": False},
        {"protein": "P1", "accession": "cl9", "name": "SF", "qstart": 1, "qend": 200,
         "evalue": 1e-10, "bitscore": 60, "superfamily": True},
    ]
    kept = dedupe(hits)
    assert [h["accession"] for h in kept] == ["cd1"]


def test_summarise_positions_and_coverage():
    hits = [
        {"protein": "P1", "accession": "cd1", "name": "A", "qstart": 1, "qend": 100,
         "evalue": 1e-30, "bitscore": 100, "superfamily": False},
        {"protein": "P1", "accession": "cd2", "name": "B", "qstart": 150, "qend": 200,
         "evalue": 1e-20, "bitscore": 80, "superfamily": False},
    ]
    rows = summarise(hits, {"P1": 300})
    assert rows[0][1] == 300 and rows[0][2] == 2
    assert "cd1:1-100" in rows[0][4]
    assert rows[0][5] == "50.3"
