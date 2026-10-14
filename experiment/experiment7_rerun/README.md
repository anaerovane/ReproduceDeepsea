# Experiment 7 — Known-variant directional check

This folder contains ref/alt inference for three reported FOXA1, GATA1, and FOXA2 mechanisms.

## Files and destinations

- `run_experiment7.py`: inference script. It uses the shared predictor names, local hg19 reference sequence files, and pretrained/epoch-53 checkpoints resolved through `model_assets.py`.
- `experiment7_raw_predictions.npz` and `experiment7_audit.json`: raw predictions, per-factor outputs, sequence windows, checkpoint hashes, and coordinate audit.
- [`RESULTS.md`](RESULTS.md): concise result table and scope.
- Code and result files are available in the matching [GitHub experiment folder](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/experiment7_rerun) and [Hugging Face experiment folder](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/experiment7_rerun).

Restore hg19 FASTA and model checkpoints using [`../DOWNLOAD.md`](../../DOWNLOAD.md).
