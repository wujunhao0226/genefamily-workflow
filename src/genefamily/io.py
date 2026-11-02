"""Small I/O helpers shared by the pipeline."""

from __future__ import annotations

import csv
import os
from typing import Dict, Iterable, List, Sequence


def read_fasta(path: str) -> Dict[str, str]:
    """Return ``{record_id: sequence}``; the first whitespace token is the id."""
    records: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_fasta(records: Iterable[tuple], path: str, width: int = 60) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def read_ids(path: str) -> List[str]:
    ids = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line.split()[0])
    return ids


def read_tsv(path: str) -> List[dict]:
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def write_tsv(path: str, header: Sequence[str], rows: Iterable[Sequence]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
