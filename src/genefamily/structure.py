"""Exon/intron structure from GFF3 annotations."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List


def parse_attributes(field: str) -> Dict[str, str]:
    attrs = {}
    for part in field.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def read_gff3(path: str) -> dict:
    genes: Dict[str, dict] = {}
    transcripts: Dict[str, dict] = {}
    exons: Dict[str, List[dict]] = defaultdict(list)

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _source, ftype, start, end, _score, strand, _phase, attr = f[:9]
            attrs = parse_attributes(attr)
            feat_id, parent = attrs.get("ID"), attrs.get("Parent")
            record = {
                "id": feat_id, "chrom": chrom, "start": int(start),
                "end": int(end), "strand": strand,
            }
            if ftype == "gene":
                genes[feat_id] = record
            elif ftype in ("mRNA", "transcript"):
                record["gene"] = parent
                transcripts[feat_id] = record
            elif ftype == "exon":
                record["transcript"] = parent
                exons[parent].append(record)
    return {"genes": genes, "transcripts": transcripts, "exons": exons}


def transcript_structure(annotation: dict) -> List[dict]:
    """Return one row per transcript with exon/intron coordinates and lengths."""
    rows: List[dict] = []
    for tid, tx in annotation["transcripts"].items():
        gene = annotation["genes"].get(tx["gene"], {})
        ex = sorted(annotation["exons"].get(tid, []), key=lambda e: e["start"])
        if not ex:
            continue
        lengths = [e["end"] - e["start"] + 1 for e in ex]
        introns = []
        for left, right in zip(ex, ex[1:]):
            introns.append({
                "start": left["end"] + 1,
                "end": right["start"] - 1,
                "length": right["start"] - left["end"] - 1,
            })
        rows.append({
            "gene_id": tx["gene"],
            "transcript_id": tid,
            "chrom": tx["chrom"],
            "strand": tx["strand"],
            "start": tx["start"],
            "end": tx["end"],
            "gene_length": gene.get("end", tx["end"]) - gene.get("start", tx["start"]) + 1,
            "n_exons": len(ex),
            "n_introns": len(introns),
            "mrna_length": sum(lengths),
            "exons": lengths,
            "introns": [i["length"] for i in introns],
        })
    return sorted(rows, key=lambda r: (r["chrom"], r["start"]))


def summary_rows(structures: List[dict]) -> List[List]:
    rows = []
    for s in structures:
        rows.append([
            s["gene_id"], s["transcript_id"], s["chrom"], s["strand"],
            s["start"], s["end"], s["gene_length"], s["n_exons"],
            s["n_introns"], s["mrna_length"],
            ",".join(map(str, s["exons"])), ",".join(map(str, s["introns"])),
        ])
    return rows


SUMMARY_HEADER = ["gene_id", "transcript_id", "chrom", "strand", "start", "end",
                  "gene_length_bp", "n_exons", "n_introns", "mrna_length_bp",
                  "exon_lengths_bp", "intron_lengths_bp"]
