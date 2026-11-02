"""Expression matrix handling: TPM normalisation and row Z-scores."""

from __future__ import annotations

from typing import List

import numpy as np

from .io import read_tsv


def read_counts(path: str, gene_col: str = "gene_id",
                length_col: str = "length") -> tuple:
    """Read a counts table; returns (genes, samples, counts matrix, lengths)."""
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"{path} is empty")
    columns = list(rows[0].keys())
    samples = [c for c in columns if c not in (gene_col, length_col)]
    genes = [r[gene_col] for r in rows]
    counts = np.array([[float(r[s]) for s in samples] for r in rows], dtype=float)
    lengths = np.array([float(r.get(length_col, 1) or 1) for r in rows], dtype=float)
    return genes, samples, counts, lengths


def tpm(counts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """TPM using the full matrix as the library (correct denominator)."""
    rpk = counts / np.maximum(lengths, 1.0)[:, None] * 1e3
    totals = rpk.sum(axis=0)
    totals[totals == 0] = 1.0
    return rpk / totals * 1e6


def zscore(matrix: np.ndarray) -> np.ndarray:
    """Row-wise Z-score with a guard for constant rows."""
    mean = matrix.mean(axis=1, keepdims=True)
    sd = matrix.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (matrix - mean) / sd


def summarise(genes: List[str], samples: List[str], tpm_matrix: np.ndarray) -> List[List]:
    rows = []
    for i, gene in enumerate(genes):
        values = tpm_matrix[i]
        best = int(np.argmax(values))
        rows.append([
            gene, f"{values.mean():.2f}", f"{values.max():.2f}", samples[best],
            int((values >= 1).sum()), len(values),
        ])
    return rows


SUMMARY_HEADER = ["gene_id", "mean_tpm", "max_tpm", "max_sample",
                  "n_samples_tpm_ge1", "n_samples"]
