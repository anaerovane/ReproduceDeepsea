# Supplementary Figure 7 feature ablation

[GitHub code, final figures, and method notes](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure7) · [Hugging Face result tables and fold caches](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure7/results)

This experiment covers the planned GRASP eQTL and GWAS cohorts. HGMD is outside
its scope, and no ClinVar proxy is used. The output consists of separate
`pretrained` and `ours` figures.

## Inputs

- Validated GRASP and GWAS variant rows, fold assignments, labels, and DeepSEA
  effect matrices from the Hugging Face Figure 3 and Experiment 9 outputs.
- Four conservation annotations (`priPhCons`, `priPhyloP`, `GerpN`, `GerpS`)
  from the completed Experiment 9 CADD GRCh37-v1.4 `inclAnno` batches.
- Only the common complete cases for all four conservation fields are scored,
  so the three feature sets are compared on the same samples. All available
  random-negative SNPs are retained (no 5% evaluation subsample).

## Recomputed models

The script fits weighted XGBoost `gblinear` classifiers using the supplied
spatial 10-fold assignments. Each fold's scaler is fit on training rows only;
training uses positives and random-negative SNPs, while evaluation uses the
out-of-fold scores against each negative group.

- `combined`: 919 DeepSEA per-track effect magnitudes plus the four CADD
  conservation annotations.
- `chromatin`: the 919 DeepSEA per-track effect magnitudes only.
- `conservation`: the four CADD conservation annotations only.

Effect magnitude follows Experiment 9:
`abs(P_ref - P_alt) * abs(logit(P_ref) - logit(P_alt))`.

## Outputs

- `suppfig7_pretrained.png` / `.pdf` (this folder)
- `suppfig7_ours.png` / `.pdf` (this folder)
- `auc_by_group.csv`: dynamically computed OOF AUCs, sample counts, and actual
  mean negative-to-positive distances.
- `run_manifest.json`: input hashes, model settings, and output hashes.
- Per-fold `.npz` predictions and `*_oof.tsv.gz` model predictions support
  auditing and resuming the run.

The final `results/` directory is published on Hugging Face. When the inputs
are absent locally, the script downloads the validated Experiment 9 cohorts
and Figure 3 effect matrices from Hugging Face; those large inputs are not
duplicated in this folder. The final PNG/PDF plots remain beside the script.

Run from the GitHub repository root with:

```bash
python -u experiment/Supplementary_Figure7/suppfig7_feature_ablation.py
```

This reproduces the planned GRASP/GWAS scope of Supplementary Figure 7 using
CADD-provided conservation annotations. Conservation annotation provenance is
not guaranteed to be identical to the original study. Results are recalculated from
local variant-level features and fold predictions, not copied from the paper.
The required cohort tables and effect matrices are fetched from the linked Hugging Face paths when they are not present locally.
