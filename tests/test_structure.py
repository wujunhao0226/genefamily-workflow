from genefamily.structure import read_gff3, transcript_structure

GFF = """##gff-version 3
chr1\tsynth\tgene\t100\t900\t.\t+\t.\tID=G1
chr1\tsynth\tmRNA\t100\t900\t.\t+\t.\tID=G1.1;Parent=G1
chr1\tsynth\texon\t100\t200\t.\t+\t.\tID=G1.1.exon1;Parent=G1.1
chr1\tsynth\texon\t400\t900\t.\t+\t.\tID=G1.1.exon2;Parent=G1.1
"""


def test_exon_and_intron_geometry(tmp_path):
    path = tmp_path / "a.gff3"
    path.write_text(GFF)
    annotation = read_gff3(str(path))
    assert set(annotation["genes"]) == {"G1"}
    rows = transcript_structure(annotation)
    assert len(rows) == 1
    row = rows[0]
    assert row["n_exons"] == 2
    assert row["n_introns"] == 1
    assert row["exons"] == [101, 501]
    assert row["introns"] == [199]          # 400 - 200 - 1
    assert row["gene_length"] == 801
