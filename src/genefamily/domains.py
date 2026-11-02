"""Conserved-domain tables: parse RPS-BLAST output, remove redundancy, summarise."""

from __future__ import annotations

import re
from typing import Dict, List

from .io import write_tsv

COLUMNS = ["qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
           "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

ACC_PREFIX = re.compile(r"^(cd|cl|COG|PRK|pfam|smart|TIGR|KOG)", re.I)


def parse_rpsblast(path: str, max_evalue: float = 1e-3) -> List[dict]:
    """Parse a 13-column RPS-BLAST tabular file (``-outfmt 6 qseqid sseqid stitle ...``)."""
    hits: List[dict] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            evalue = float(f[11])
            if evalue > max_evalue:
                continue
            acc, name = _identity(f[1], f[2])
            hits.append({
                "protein": f[0],
                "accession": acc,
                "name": name,
                "qstart": int(f[7]),
                "qend": int(f[8]),
                "evalue": evalue,
                "bitscore": float(f[12]),
                "superfamily": acc.lower().startswith("cl"),
            })
    return hits


def _identity(sseqid: str, stitle: str):
    """RPS-BLAST reports a numeric UID; the real model id/name live in ``stitle``."""
    title = stitle.strip()
    if ACC_PREFIX.match(title):
        parts = [p.strip() for p in title.split(",")]
        return parts[0], (parts[1] if len(parts) > 1 else parts[0])
    return sseqid.split("|")[-1].strip(), title


def _overlap(a: dict, b: dict) -> int:
    return max(0, min(a["qend"], b["qend"]) - max(a["qstart"], b["qstart"]) + 1)


def dedupe(hits: List[dict]) -> List[dict]:
    """Keep the best specific hit per overlapping region; drop redundant superfamily hits."""
    by_protein: Dict[str, List[dict]] = {}
    for hit in hits:
        by_protein.setdefault(hit["protein"], []).append(hit)

    kept_all: List[dict] = []
    for protein, group in by_protein.items():
        specific = sorted([h for h in group if not h["superfamily"]],
                          key=lambda h: (h["qstart"], -h["bitscore"]))
        superfamily = [h for h in group if h["superfamily"]]

        kept: List[dict] = []
        for hit in specific:
            clash = None
            for other in kept:
                shorter = min(hit["qend"] - hit["qstart"] + 1,
                              other["qend"] - other["qstart"] + 1)
                if _overlap(hit, other) >= 0.5 * shorter:
                    clash = other
                    break
            if clash is None:
                kept.append(hit)
            elif (hit["evalue"], -hit["bitscore"]) < (clash["evalue"], -clash["bitscore"]):
                kept[kept.index(clash)] = hit

        for hit in superfamily:
            if any(_overlap(hit, k) >= 0.5 * min(hit["qend"] - hit["qstart"] + 1,
                                                k["qend"] - k["qstart"] + 1)
                   for k in kept):
                continue
            kept.append(hit)

        kept_all.extend(sorted(kept, key=lambda h: h["qstart"]))
    return kept_all


def summarise(hits: List[dict], lengths: Dict[str, int]) -> List[List]:
    """One row per protein: number of domains, types and their positions."""
    by_protein: Dict[str, List[dict]] = {}
    for hit in hits:
        by_protein.setdefault(hit["protein"], []).append(hit)

    rows = []
    for protein, group in sorted(by_protein.items()):
        group = sorted(group, key=lambda h: h["qstart"])
        types = ";".join(sorted({f'{h["accession"]}:{h["name"]}' for h in group}))
        positions = ";".join(f'{h["accession"]}:{h["qstart"]}-{h["qend"]}' for h in group)
        covered = sum(h["qend"] - h["qstart"] + 1 for h in group)
        length = lengths.get(protein, 0)
        rows.append([
            protein, length, len(group), types, positions,
            f"{100.0 * covered / length:.1f}" if length else "0.0",
        ])
    return rows


def write_tables(hits: List[dict], outdir: str, prefix: str = "domains") -> dict:
    hits_tsv = f"{outdir}/{prefix}_hits.tsv"
    write_tsv(hits_tsv,
              ["protein", "accession", "name", "start", "end", "evalue", "bitscore", "class"],
              [[h["protein"], h["accession"], h["name"], h["qstart"], h["qend"],
                f'{h["evalue"]:.2e}', f'{h["bitscore"]:.1f}',
                "superfamily" if h["superfamily"] else "specific"] for h in hits])
    best = dedupe(hits)
    best_tsv = f"{outdir}/{prefix}_hits_best.tsv"
    write_tsv(best_tsv,
              ["protein", "accession", "name", "start", "end", "evalue", "bitscore", "class"],
              [[h["protein"], h["accession"], h["name"], h["qstart"], h["qend"],
                f'{h["evalue"]:.2e}', f'{h["bitscore"]:.1f}',
                "superfamily" if h["superfamily"] else "specific"] for h in best])
    return {"hits": hits_tsv, "best": best_tsv}
