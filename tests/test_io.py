from genefamily.io import read_fasta, read_ids, write_fasta


def test_fasta_roundtrip(tmp_path):
    path = tmp_path / "x.faa"
    write_fasta([("a", "MKV"), ("b", "MKKLL")], str(path), width=2)
    records = read_fasta(str(path))
    assert records == {"a": "MKV", "b": "MKKLL"}
    text = path.read_text()
    assert text.startswith(">a\nMK\nV\n")


def test_read_ids_skips_comments(tmp_path):
    path = tmp_path / "ids.txt"
    path.write_text("# header\nG1\nG2 extra\n\n")
    assert read_ids(str(path)) == ["G1", "G2"]
