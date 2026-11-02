"""Publication-style figures: heatmap, circular tree and gene models."""

from __future__ import annotations

import math
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from .tree import all_nodes, leaves, parse_newick  # noqa: E402

CMAP = LinearSegmentedColormap.from_list(
    "RdBu_trunc",
    ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
     "#fddbc7", "#f4a582", "#d6604d", "#b2182b"],
)
PALETTE = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
           "#76B7B2", "#EDC948", "#9C755F"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def save(fig, path: str, formats: List[str], dpi: int = 300) -> None:
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"  wrote {path}.{fmt}")
    plt.close(fig)


def heatmap(matrix: np.ndarray, row_labels: List[str], col_labels: List[str],
            path: str, title: str, formats: List[str], dpi: int = 300,
            vmax: float = 2.0) -> None:
    n_rows, n_cols = matrix.shape
    fig, ax = plt.subplots(figsize=(1.0 + 0.45 * n_cols, 0.25 * n_rows + 2.0))
    fig.subplots_adjust(left=0.28, right=0.88, top=0.88, bottom=0.18)
    im = ax.pcolormesh(np.arange(n_cols + 1), np.arange(n_rows + 1), matrix,
                       cmap=CMAP, vmin=-vmax, vmax=vmax,
                       edgecolors="#4a5568", linewidth=0.3, shading="flat")
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(row_labels, fontsize=6, family="DejaVu Sans Mono")
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Row Z-score of log2(TPM+1)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    save(fig, path, formats, dpi)


def circular_tree(newick: str, groups: Dict[str, str], path: str, title: str,
                  formats: List[str], dpi: int = 300) -> None:
    root = parse_newick(newick)
    tips = leaves(root)
    n = len(tips)

    angles, dist = {}, {}
    for i, tip in enumerate(tips):
        angles[id(tip)] = 2 * math.pi * (i + 0.5) / n
    dist[id(root)] = 0.0
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            dist[id(child)] = dist[id(node)] + max(child.length, 0.0)
            stack.append(child)
    dmax = max(dist.values()) or 1.0

    def angle(node):
        if not node.children:
            return angles[id(node)]
        value = sum(angle(c) for c in node.children) / len(node.children)
        angles[id(node)] = value
        return value

    angle(root)
    group_names = sorted(set(groups.get(t.name, "NA") for t in tips))
    colors = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(group_names)}

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.38, 1.38)
    ax.set_ylim(-1.38, 1.38)

    def radius(node):
        if not node.children:
            return 1.0
        return 0.05 + 0.90 * (dist[id(node)] / dmax) ** 0.5

    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            a = angles[id(child)]
            r0, r1 = radius(node), radius(child)
            color = colors.get(groups.get(child.name, "NA"), "#9aa0a6") if not child.children \
                else "#b8bec5"
            ax.plot([r0 * math.cos(a), r1 * math.cos(a)],
                    [r0 * math.sin(a), r1 * math.sin(a)],
                    color=color, lw=0.9 if not child.children else 0.7)
            if child.children:
                span = 2 * math.pi / n * (len(leaves(child)) - 1)
                theta = np.linspace(a - span / 2, a + span / 2, 24)
                ax.plot(r1 * np.cos(theta), r1 * np.sin(theta), color="#b8bec5", lw=0.7)
                stack.append(child)

    for node in all_nodes(root):
        if node.children and node.support is not None:
            a = angles[id(node)]
            r = radius(node)
            ax.scatter([r * math.cos(a)], [r * math.sin(a)], s=6 if node.support >= 95 else 4,
                       c="#1a1a1a" if node.support >= 95 else "#8a8f94", zorder=5)

    for tip in tips:
        a = angles[id(tip)]
        deg = math.degrees(a) % 360
        if 90 <= deg <= 270:
            rot, ha = deg + 180, "right"
        else:
            rot, ha = deg, "left"
        ax.text(1.03 * math.cos(a), 1.03 * math.sin(a), tip.name,
                rotation=rot, rotation_mode="anchor", ha=ha, va="center",
                fontsize=6, color=colors.get(groups.get(tip.name, "NA"), "#4a4a4a"))

    handles = [plt.Line2D([0], [0], marker="o", color="none",
                          markerfacecolor=colors[g], markersize=5, label=g)
               for g in group_names]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 0.0),
              frameon=False, fontsize=7, title="Group", title_fontsize=8)
    ax.set_title(title, fontsize=11, fontweight="bold")
    save(fig, path, formats, dpi)


def gene_structure(structures: List[dict], path: str, title: str,
                   formats: List[str], dpi: int = 300) -> None:
    n = len(structures)
    fig, ax = plt.subplots(figsize=(10, 0.32 * n + 1.6))
    fig.subplots_adjust(left=0.28, right=0.96, top=0.90, bottom=0.12)
    max_len = max(s["gene_length"] for s in structures) or 1

    for i, s in enumerate(structures):
        y = n - 1 - i
        ax.plot([0, s["gene_length"]], [y, y], color="#4a5568", lw=0.8, zorder=1)
        offset = 0
        for length in s["exons"]:
            ax.add_patch(Rectangle((offset, y - 0.22), max(length, 12), 0.44,
                                   facecolor="#2b6cb0", edgecolor="none", zorder=3))
            offset += length + (s["introns"][0] if s["introns"] else 0)
        ax.text(-0.01 * max_len, y, s["gene_id"], transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6.5, family="DejaVu Sans Mono")
        ax.text(1.01, y, f'{s["n_exons"]} exons', transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=6, color="#4a5568")

    ax.set_xlim(0, max_len)
    ax.set_ylim(-0.8, n + 0.4)
    ax.set_yticks([])
    ax.set_xlabel("Gene position (bp, schematic)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    save(fig, path, formats, dpi)
