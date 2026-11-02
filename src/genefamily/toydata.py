"""Deterministic synthetic dataset used by the tests and by CI.

The sequences are random strings with a planted motif; they are **not** biological
data and carry no information from any real genome or experiment.
"""

from __future__ import annotations

import csv
import os
import random
from typing import List

AA = "ACDEFGHIKLMNPQRSTVWY"
MOTIF = ("MEEVCGCEGPALGQRLVAALRARGIDVSELEQLRAELGRLRRRLQAAGEV")
SPECIES = ["spA", "spB", "spC"]
GROUPS = ["alpha", "beta", "gamma"]


def _random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(AA) for _ in range(length))


def _mutate(rng: random.Random, seq: str, rate: float = 0.08) -> str:
    return "".join(rng.choice(AA) if rng.random() < rate else c for c in seq)


def build(outdir: str = "data/toy", seed: int = 42, per_species: int = 20) -> dict:
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)

    records, groups, hit_rows, gff_lines = [], [], [], ["##gff-version 3"]
    motif_ids: List[str] = []

    for sp in SPECIES:
        for i in range(1, per_species + 1):
            gid = f"TOY|{sp}|G{i:04d}"
            has_motif = i <= 8
            core = _random_seq(rng, rng.randint(60, 100))
            if has_motif:
                pos = rng.randint(20, 60)
                core = core[:pos] + _mutate(rng, MOTIF) + core[pos:]
                motif_ids.append(gid)
                hit_rows.append((gid, f"cd{20000 + i:05d}"))
            records.append((gid, core))
            groups.append((gid, sp, GROUPS[i % 3]))

            # simple three-exon model in GFF3
            start = 1000 + i * 5000
            exon_len = max(120, len(core) // 3)
            exons = [(start, start + exon_len - 1),
                     (start + exon_len + 300, start + 2 * exon_len + 299),
                     (start + 2 * exon_len + 700, start + 3 * exon_len + 699)]
            gff_lines.append(f"{sp}\tsynth\tgene\t{start}\t{exons[-1][1]}\t.\t+\t.\tID={gid}")
            gff_lines.append(f"{sp}\tsynth\tmRNA\t{start}\t{exons[-1][1]}\t.\t+\t.\tID={gid}.1;Parent={gid}")
            for k, (s, e) in enumerate(exons, 1):
                gff_lines.append(
                    f"{sp}\tsynth\texon\t{s}\t{e}\t.\t+\t.\tID={gid}.1.exon{k};Parent={gid}.1")

    with open(os.path.join(outdir, "proteome.faa"), "w", encoding="utf-8") as fh:
        for gid, seq in records:
            fh.write(f">{gid} synthetic=true\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    with open(os.path.join(outdir, "seeds.faa"), "w", encoding="utf-8") as fh:
        for k in range(8):
            fh.write(f">seed{k}\n{_mutate(rng, MOTIF, 0.03)}\n")

    with open(os.path.join(outdir, "annotation.gff3"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(gff_lines) + "\n")

    with open(os.path.join(outdir, "groups.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "species", "group"])
        w.writerows(groups)

    # synthetic RPS-BLAST-like output (13 tab-separated columns)
    with open(os.path.join(outdir, "rpsblast_hits.tsv"), "w", encoding="utf-8") as fh:
        for gid, acc in hit_rows:
            fh.write("\t".join([
                gid, f"gnl|CDD|{acc[2:]}", f"{acc}, CYCLIN_TOY, synthetic model",
                "31.0", "120", "10", "2", "5", "145", "1", "120", "1e-30", "120.0",
            ]) + "\n")
            fh.write("\t".join([
                gid, f"gnl|CDD|{int(acc[2:]) + 1}", f"cd{int(acc[2:]) + 1:05d}, CYCLIN_TOY_C, synthetic model",
                "28.0", "90", "8", "1", "160", "250", "1", "90", "2e-18", "80.0",
            ]) + "\n")

    # synthetic expression matrix: motif genes are group-specific
    samples = [f"{g}_{r}" for g in GROUPS for r in range(1, 4)]
    with open(os.path.join(outdir, "counts.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "length"] + samples)
        for gid, seq in records:
            base = rng.randint(10, 80)
            row = []
            for sample in samples:
                group = sample.split("_")[0]
                if gid in motif_ids and group == GROUPS[hash(gid) % 3]:
                    value = base * rng.randint(20, 60)
                else:
                    value = base * rng.randint(1, 3)
                row.append(value)
            w.writerow([gid, len(seq) * 3] + row)

    # synthetic tree so that plotting can be exercised without IQ-TREE
    tips = [g for g, _ in records][:12]
    tree = _balanced_newick(tips)
    with open(os.path.join(outdir, "toy_tree.nwk"), "w", encoding="utf-8") as fh:
        fh.write(tree + "\n")

    return {"outdir": outdir, "n_sequences": len(records), "n_motif": len(motif_ids)}


def _balanced_newick(names: List[str]) -> str:
    if len(names) == 1:
        return f"{names[0]}:0.1"

    def build(chunk):
        if len(chunk) == 1:
            return f"{chunk[0]}:0.1"
        mid = len(chunk) // 2
        return f"({build(chunk[:mid])},{build(chunk[mid:])})90:0.05"

    return build(names) + ";"
