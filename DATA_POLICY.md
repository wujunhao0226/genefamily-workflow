# Data policy

This repository ships **no experimental data**.

- Everything under `data/toy/` is generated at runtime by `scripts/make_toy_data.py`
  from a fixed random seed. Sequences are artificial, identifiers are artificial
  (`TOY...`), and the "species" are placeholders. They carry no biological meaning.
- Real input data must be placed under `data/raw/`, which is git-ignored
  (see `.gitignore`). Nothing in that directory can be committed by accident.
- Generated outputs go to `results/`, also git-ignored.
- The CI workflow only ever executes the synthetic toy dataset, so no private
  sequence, sample name, coordinate or expression value can appear in build logs.

If you fork this repository for your own project, keep the same layout:

```
data/raw/      <- your data, never committed
data/toy/      <- synthetic, regenerated on demand
results/       <- pipeline output, never committed
```
