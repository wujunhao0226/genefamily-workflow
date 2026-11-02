"""End-to-end run on the synthetic dataset (no external binaries required)."""

import os

from genefamily import cli


def test_toy_dataset_to_figures(tmp_path):
    data = tmp_path / "data"
    results = tmp_path / "results"
    assert cli.main(["toy", "--outdir", str(data), "--per-species", "12"]) == 0

    assert cli.main(["structure", "--gff", str(data / "annotation.gff3"),
                     "--outdir", str(results)]) == 0
    assert cli.main(["expression", "--counts", str(data / "counts.tsv"),
                     "--outdir", str(results)]) == 0
    assert cli.main(["domains", "--hits", str(data / "rpsblast_hits.tsv"),
                     "--fasta", str(data / "proteome.faa"),
                     "--outdir", str(results)]) == 0
    assert cli.main(["figures", "--results", str(results),
                     "--tree", str(data / "toy_tree.nwk"),
                     "--groups", str(data / "groups.tsv"),
                     "--gff", str(data / "annotation.gff3"),
                     "--outdir", str(results / "figures"),
                     "--formats", "png", "--dpi", "120"]) == 0

    for name in ("gene_structure.tsv", "tpm.tsv", "zscore.tsv",
                 "expression_summary.tsv", "domain_summary.tsv"):
        assert (results / name).is_file(), name

    figures = os.listdir(results / "figures")
    assert any(f.startswith("circular_tree") for f in figures)
    assert any(f.startswith("expression_heatmap") for f in figures)
    assert any(f.startswith("gene_models") for f in figures)
