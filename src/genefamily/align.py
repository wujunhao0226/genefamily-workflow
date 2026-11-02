"""Multiple sequence alignment wrappers (Clustal Omega or MAFFT)."""

from __future__ import annotations

import os
import subprocess

from .hmmer import ToolMissingError, require


def run_alignment(fasta: str, out_fasta: str, tool: str = "auto",
                  threads: int = 4) -> str:
    """Align ``fasta`` and write ``out_fasta``.

    ``tool`` accepts ``clustalo``, ``mafft`` or ``auto`` (first one found).
    """
    import shutil

    if tool == "auto":
        for candidate in ("clustalo", "mafft"):
            if shutil.which(candidate):
                tool = candidate
                break
        else:
            raise ToolMissingError(
                "Neither 'clustalo' nor 'mafft' is on PATH. "
                "Install one (`micromamba install -c bioconda clustalo`) or use --skip-external."
            )

    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)), exist_ok=True)
    require(tool)

    if tool == "clustalo":
        cmd = [tool, "-i", fasta, "-o", out_fasta, "--outfmt=fasta",
               "--seqtype=protein", "--force", "--threads=%d" % threads]
    elif tool == "mafft":
        cmd = [tool, "--auto", "--thread", str(threads), fasta]
        subprocess.run(cmd, check=True, stdout=open(out_fasta, "w", encoding="utf-8"))
        return out_fasta
    else:
        raise ValueError(f"unsupported aligner: {tool}")

    subprocess.run(cmd, check=True)
    return out_fasta


def alignment_stats(fasta: str) -> dict:
    """Basic sanity metrics: number of sequences, columns, gap fraction."""
    from .io import read_fasta

    records = read_fasta(fasta)
    if not records:
        return {"n_sequences": 0, "n_columns": 0, "gap_fraction": 0.0}
    lengths = {len(s) for s in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"alignment has unequal lengths: {sorted(lengths)[:5]}")
    ncol = lengths.pop()
    gaps = sum(s.count("-") for s in records.values())
    return {
        "n_sequences": len(records),
        "n_columns": ncol,
        "gap_fraction": gaps / (len(records) * ncol) if ncol else 0.0,
    }
