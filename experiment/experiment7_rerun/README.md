# Experiment 7 — Known-variant directional check

This folder contains a three-locus ref/alt inference check for the reported FOXA1, GATA1, and FOXA2 mechanisms. It is a small directional spot-check, not a full cohort reproduction.

## Files and destinations

- `run_experiment7.py`: inference script. It uses the shared predictor names, local hg19 reference sequence files, and pretrained/epoch-53 checkpoints resolved through `model_assets.py`.
- `experiment7_raw_predictions.npz` and `experiment7_audit.json`: raw predictions, per-factor outputs, sequence windows, checkpoint hashes, and coordinate audit.
- [`RESULTS.md`](RESULTS.md): concise result table and scope.
- Code and these small result files are published only in the matching [GitHub experiment folder](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/experiment7_rerun). No separate Hugging Face artifact is needed for these small files.

Restore hg19 FASTA and model checkpoints using [`../DOWNLOAD.md`](../../DOWNLOAD.md); no reference genome or model copy is stored in this experiment folder.
