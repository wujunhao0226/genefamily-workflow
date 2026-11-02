"""HMMER discovery step: run ``hmmsearch`` and parse its tabular output."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from typing import Dict, List

from .io import read_fasta, write_fasta

TOOL = "hmmsearch"


class ToolMissingError(RuntimeError):
    """Raised when an external binary is not on PATH."""


def require(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        raise ToolMissingError(
            f"'{tool}' not found on PATH. Install it (conda: `micromamba install -c bioconda {tool}`) "
            f"or skip this step with --skip-external."
        )
    return path


def run_hmmbuild(seeds: str, hmm: str) -> str:
    """Build a profile HMM from seed sequences."""
    require("hmmbuild")
    os.makedirs(os.path.dirname(os.path.abspath(hmm)), exist_ok=True)
    subprocess.run(["hmmbuild", "--amino", "-o", os.devnull, hmm, seeds], check=True)
    return hmm


def run_hmmsearch(hmm: str, proteome: str, tblout: str, threads: int = 4,
                  cpu: bool = True) -> str:
    """Run hmmsearch and return the path to the tabular output."""
    require(TOOL)
    os.makedirs(os.path.dirname(os.path.abspath(tblout)), exist_ok=True)
    cmd = [TOOL, "--tblout", tblout, "--noali", "-o", os.devnull]
    if cpu:
        cmd += ["--cpu", str(threads)]
    cmd += [hmm, proteome]
    subprocess.run(cmd, check=True)
    return tblout


def parse_tblout(path: str) -> List[dict]:
    """Parse a HMMER3 ``--tblout`` file into a list of hit dictionaries."""
    hits: List[dict] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 6:
                continue
            hits.append({
                "target": f[0],
                "target_acc": "-" if f[1] == "-" else f[1],
                "query": f[2],
                "query_acc": f[3],
                "evalue": float(f[4]),
                "score": float(f[5]),
                "description": " ".join(f[18:]) if len(f) > 18 else "",
            })
    return hits


def select_hits(hits: List[dict], max_evalue: float = 1e-5,
                best_per_gene: bool = True) -> List[dict]:
    """Filter by E-value and optionally collapse isoforms to the best hit per gene."""
    kept = [h for h in hits if h["evalue"] <= max_evalue]
    kept.sort(key=lambda h: (h["target"], h["evalue"]))
    if not best_per_gene:
        return kept
    best: Dict[str, dict] = {}
    for h in kept:
        gene = re.sub(r"\.\d+$", "", h["target"])
        if gene not in best or h["evalue"] < best[gene]["evalue"]:
            best[gene] = h
    return sorted(best.values(), key=lambda h: h["evalue"])


def extract_sequences(proteome: str, hits: List[dict], out_fasta: str) -> str:
    records = read_fasta(proteome)
    wanted = [h["target"] for h in hits]
    missing = [w for w in wanted if w not in records]
    if missing:
        raise KeyError(f"{len(missing)} hit(s) absent from proteome, e.g. {missing[:3]}")
    write_fasta(((w, records[w]) for w in wanted), out_fasta)
    return out_fasta
