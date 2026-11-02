"""IQ-TREE wrapper and a minimal Newick parser used by the plotting code."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from typing import List, Optional

from .hmmer import ToolMissingError


def _iqtree_binary() -> str:
    for name in ("iqtree2", "iqtree"):
        if shutil.which(name):
            return name
    raise ToolMissingError(
        "IQ-TREE not found on PATH. Install it (`micromamba install -c bioconda iqtree`) "
        "or run with --skip-external."
    )


def run_iqtree(alignment: str, prefix: str, model: str = "MFP",
               ufboot: int = 1000, threads: str = "AUTO",
               extra: Optional[List[str]] = None) -> dict:
    """Run IQ-TREE and return the paths of the produced files."""
    binary = _iqtree_binary()
    os.makedirs(os.path.dirname(os.path.abspath(prefix)) or ".", exist_ok=True)
    cmd = [binary, "-s", alignment, "-m", model, "-B", str(ufboot),
           "-T", str(threads), "--prefix", prefix, "--redo", "--quiet"]
    if extra:
        cmd += extra
    subprocess.run(cmd, check=True)
    return {
        "treefile": f"{prefix}.treefile",
        "report": f"{prefix}.iqtree",
        "log": f"{prefix}.log",
    }


class Node:
    __slots__ = ("name", "length", "support", "children", "parent")

    def __init__(self) -> None:
        self.name = ""
        self.length = 0.0
        self.support: Optional[float] = None
        self.children: List["Node"] = []
        self.parent: Optional["Node"] = None


def parse_newick(text: str) -> Node:
    """Parse a Newick string with branch lengths and numeric support values."""
    text = text.strip().rstrip(";")
    pos = 0

    def parse() -> Node:
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                child = parse()
                child.parent = node
                node.children.append(child)
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"malformed Newick at position {pos}")
        m = re.match(r"([^,():;]+)", text[pos:])
        label = ""
        if m:
            label = m.group(1).strip()
            pos += len(m.group(1))
        if text[pos:pos + 1] == ":":
            pos += 1
            m = re.match(r"([0-9.eE+-]+)", text[pos:])
            node.length = float(m.group(1))
            pos += len(m.group(1))
        if node.children and re.fullmatch(r"[0-9.]+", label or "x"):
            node.support = float(label)
        else:
            node.name = label
        return node

    root = parse()
    root.length = 0.0
    return root


def leaves(root: Node) -> List[Node]:
    if not root.children:
        return [root]
    out: List[Node] = []
    for child in root.children:
        out.extend(leaves(child))
    return out


def all_nodes(root: Node) -> List[Node]:
    out, stack = [], [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(node.children)
    return out
