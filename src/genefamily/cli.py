"""Command line interface for the gene-family workflow."""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import yaml

from . import align, domains, expression, hmmer, plotting, structure, toydata
from .io import read_fasta, read_ids, read_tsv, write_tsv
from .structure import read_gff3, transcript_structure
from .tree import leaves, parse_newick, run_iqtree


def _log(message: str) -> None:
    print(f"[genefamily] {message}", flush=True)


def _ensure_dirs(*paths: str) -> None:
    for path in paths:
        os.makedirs(path, exist_ok=True)


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------

def cmd_toy(args) -> int:
    info = toydata.build(args.outdir, seed=args.seed, per_species=args.per_species)
    _log(f"synthetic dataset: {info['n_sequences']} sequences, {info['n_motif']} with motif")
    return 0


def cmd_discover(args) -> int:
    _ensure_dirs(args.outdir)
    hmm = args.hmm
    if not os.path.isfile(hmm):
        _log(f"building HMM from {args.seeds} -> {hmm}")
        hmmer.run_hmmbuild(args.seeds, hmm)
    tblout = os.path.join(args.outdir, "hmmsearch.tblout")
    _log("running hmmsearch")
    hmmer.run_hmmsearch(hmm, args.proteome, tblout, threads=args.threads)
    hits = hmmer.parse_tblout(tblout)
    selected = hmmer.select_hits(hits, max_evalue=args.evalue,
                                 best_per_gene=not args.keep_isoforms)
    ids = [h["target"] for h in selected]
    write_tsv(os.path.join(args.outdir, "hits.tsv"),
              ["target", "query", "evalue", "score"],
              [[h["target"], h["query"], f'{h["evalue"]:.2e}', f'{h["score"]:.1f}']
               for h in selected])
    with open(os.path.join(args.outdir, "ids.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(ids) + "\n")
    fasta = os.path.join(args.outdir, "hits.faa")
    hmmer.extract_sequences(args.proteome, selected, fasta)
    _log(f"{len(hits)} raw hits -> {len(selected)} selected -> {fasta}")
    return 0


def cmd_align(args) -> int:
    _ensure_dirs(os.path.dirname(os.path.abspath(args.out)))
    align.run_alignment(args.fasta, args.out, tool=args.tool, threads=args.threads)
    stats = align.alignment_stats(args.out)
    _log(f"alignment: {stats['n_sequences']} sequences x {stats['n_columns']} columns, "
         f"gap fraction {stats['gap_fraction']:.3f}")
    return 0


def cmd_tree(args) -> int:
    _ensure_dirs(os.path.dirname(os.path.abspath(args.prefix)))
    out = run_iqtree(args.alignment, args.prefix, model=args.model,
                     ufboot=args.ufboot, threads=args.threads)
    _log(f"tree written: {out['treefile']}")
    return 0


def cmd_domains(args) -> int:
    _ensure_dirs(args.outdir)
    hits = domains.parse_rpsblast(args.hits, max_evalue=args.evalue)
    best = domains.dedupe(hits)
    lengths = {k: len(v) for k, v in read_fasta(args.fasta).items()}
    tables = domains.write_tables(hits, args.outdir)
    rows = domains.summarise(best, lengths)
    write_tsv(os.path.join(args.outdir, "domain_summary.tsv"),
              ["protein", "length", "n_domains", "domain_types", "positions", "coverage_pct"],
              rows)
    _log(f"{len(hits)} domain hits -> {len(best)} non-redundant; tables in {tables['hits']}")
    return 0


def cmd_structure(args) -> int:
    _ensure_dirs(args.outdir)
    annotation = structure.read_gff3(args.gff)
    transcripts = structure.transcript_structure(annotation)
    write_tsv(os.path.join(args.outdir, "gene_structure.tsv"),
              structure.SUMMARY_HEADER, structure.summary_rows(transcripts))
    _log(f"{len(transcripts)} transcripts from {len(annotation['genes'])} genes")
    return 0


def cmd_expression(args) -> int:
    _ensure_dirs(args.outdir)
    genes, samples, counts, lengths = expression.read_counts(args.counts)
    if args.genes:
        wanted = set(read_ids(args.genes))
        keep = [i for i, g in enumerate(genes) if g in wanted]
        genes = [genes[i] for i in keep]
        counts = counts[keep]
        lengths = lengths[keep]
    tpm = expression.tpm(counts, lengths)
    z = expression.zscore(np.log2(tpm + 1))
    write_tsv(os.path.join(args.outdir, "tpm.tsv"), ["gene_id"] + samples,
              [[g] + [f"{v:.3f}" for v in row] for g, row in zip(genes, tpm)])
    write_tsv(os.path.join(args.outdir, "zscore.tsv"), ["gene_id"] + samples,
              [[g] + [f"{v:.3f}" for v in row] for g, row in zip(genes, z)])
    write_tsv(os.path.join(args.outdir, "expression_summary.tsv"),
              expression.SUMMARY_HEADER,
              expression.summarise(genes, samples, tpm))
    _log(f"{len(genes)} genes x {len(samples)} samples -> TPM/Z-score")
    return 0


def cmd_figures(args) -> int:
    _ensure_dirs(args.outdir)
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]

    tree_path = args.tree
    if not os.path.isfile(tree_path):
        candidates = sorted(glob.glob(os.path.join(args.results, "*.treefile")))
        if not candidates:
            raise FileNotFoundError(f"no treefile found (looked for {tree_path})")
        tree_path = candidates[0]
    groups = {}
    if args.groups and os.path.isfile(args.groups):
        for row in read_tsv(args.groups):
            groups[row["gene_id"]] = row.get("group", row.get("species", "NA"))
    newick = open(tree_path, encoding="utf-8").read()
    tips = [t.name for t in leaves(parse_newick(newick))]
    if not groups:
        groups = {t: t.split("|")[1] if "|" in t else "NA" for t in tips}
    plotting.circular_tree(newick, groups, os.path.join(args.outdir, "circular_tree"),
                           "Synthetic gene-family phylogeny", formats, dpi=args.dpi)

    heat = os.path.join(args.results, "zscore.tsv")
    if os.path.isfile(heat):
        rows = read_tsv(heat)
        samples = [c for c in rows[0] if c != "gene_id"]
        matrix = np.array([[float(r[s]) for s in samples] for r in rows])
        plotting.heatmap(matrix, [r["gene_id"] for r in rows], samples,
                         os.path.join(args.outdir, "expression_heatmap"),
                         "Synthetic expression (row Z-score)", formats, dpi=args.dpi)

    struct_table = os.path.join(args.results, "gene_structure.tsv")
    if os.path.isfile(struct_table) and args.gff and os.path.isfile(args.gff):
        models = transcript_structure(read_gff3(args.gff))[:12]
        plotting.gene_structure(models, os.path.join(args.outdir, "gene_models"),
                                "Synthetic gene models", formats, dpi=args.dpi)
    return 0


def cmd_run(args) -> int:
    with open(args.config, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    results = cfg.get("results", "results")
    _ensure_dirs(results)
    skip = args.skip_external

    seeds = cfg["seeds"]
    hmm = cfg.get("hmm", os.path.join(results, "profile.hmm"))
    if not skip:
        cmd_discover(argparse.Namespace(outdir=os.path.join(results, "discover"),
                                        hmm=hmm, seeds=seeds, proteome=cfg["proteome"],
                                        evalue=cfg.get("evalue", 1e-5), threads=args.threads,
                                        keep_isoforms=False))
        cmd_align(argparse.Namespace(fasta=os.path.join(results, "discover", "hits.faa"),
                                     out=os.path.join(results, "alignment.faa"),
                                     tool=cfg.get("aligner", "auto"), threads=args.threads))
        cmd_tree(argparse.Namespace(alignment=os.path.join(results, "alignment.faa"),
                                    prefix=os.path.join(results, "tree"),
                                    model=cfg.get("model", "MFP"),
                                    ufboot=cfg.get("ufboot", 1000),
                                    threads=args.threads))
        tree_file = os.path.join(results, "tree.treefile")
    else:
        _log("skipping external steps (--skip-external)")
        tree_file = cfg.get("tree", "data/toy/toy_tree.nwk")

    if cfg.get("gff"):
        cmd_structure(argparse.Namespace(gff=cfg["gff"], outdir=results))
    if cfg.get("counts"):
        cmd_expression(argparse.Namespace(counts=cfg["counts"], genes=None,
                                          outdir=results))
    if cfg.get("domains"):
        fasta = os.path.join(results, "discover", "hits.faa")
        if not os.path.isfile(fasta):
            fasta = cfg["proteome"]
        cmd_domains(argparse.Namespace(hits=cfg["domains"], fasta=fasta,
                                       outdir=results, evalue=1e-3))

    cmd_figures(argparse.Namespace(results=results, tree=tree_file,
                                   groups=cfg.get("groups"), gff=cfg.get("gff"),
                                   outdir=os.path.join(results, "figures"),
                                   formats=cfg.get("formats", "png"),
                                   dpi=cfg.get("dpi", 300)))

    manifest = []
    for path in sorted(glob.glob(os.path.join(results, "**", "*"), recursive=True)):
        if os.path.isfile(path):
            manifest.append([os.path.relpath(path, results), os.path.getsize(path)])
    write_tsv(os.path.join(results, "manifest.tsv"), ["file", "size_bytes"], manifest)
    _log(f"pipeline finished, {len(manifest)} files under {results}")
    return 0


# --------------------------------------------------------------------------
# argument parsing
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genefamily",
        description="Reproducible gene-family workflow (HMMER -> MSA -> IQ-TREE -> figures).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("toy", help="generate the synthetic dataset")
    p.add_argument("--outdir", default="data/toy")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--per-species", type=int, default=20)
    p.set_defaults(func=cmd_toy)

    p = sub.add_parser("discover", help="HMM search and hit extraction")
    p.add_argument("--hmm", required=True)
    p.add_argument("--seeds")
    p.add_argument("--proteome", required=True)
    p.add_argument("--outdir", default="results/discover")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-isoforms", action="store_true")
    p.set_defaults(func=cmd_discover)

    p = sub.add_parser("align", help="multiple sequence alignment")
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tool", default="auto", choices=["auto", "clustalo", "mafft"])
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("tree", help="IQ-TREE phylogeny")
    p.add_argument("--alignment", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--model", default="MFP")
    p.add_argument("--ufboot", type=int, default=1000)
    p.add_argument("--threads", type=str, default="AUTO")
    p.set_defaults(func=cmd_tree)

    p = sub.add_parser("domains", help="domain tables from RPS-BLAST output")
    p.add_argument("--hits", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--outdir", default="results")
    p.add_argument("--evalue", type=float, default=1e-3)
    p.set_defaults(func=cmd_domains)

    p = sub.add_parser("structure", help="exon/intron tables from GFF3")
    p.add_argument("--gff", required=True)
    p.add_argument("--outdir", default="results")
    p.set_defaults(func=cmd_structure)

    p = sub.add_parser("expression", help="TPM and Z-score tables")
    p.add_argument("--counts", required=True)
    p.add_argument("--genes")
    p.add_argument("--outdir", default="results")
    p.set_defaults(func=cmd_expression)

    p = sub.add_parser("figures", help="render figures")
    p.add_argument("--results", default="results")
    p.add_argument("--tree", default="results/tree.treefile")
    p.add_argument("--groups")
    p.add_argument("--gff")
    p.add_argument("--outdir", default="results/figures")
    p.add_argument("--formats", default="png")
    p.add_argument("--dpi", type=int, default=300)
    p.set_defaults(func=cmd_figures)

    p = sub.add_parser("run", help="run the whole pipeline from a YAML config")
    p.add_argument("--config", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--skip-external", action="store_true")
    p.set_defaults(func=cmd_run)

    return parser


def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except hmmer.ToolMissingError as exc:
        print(f"[genefamily] missing tool: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
