# Experiment 9 — Unsupervised functional significance

The completed analysis uses GRASP eQTL and GWAS Catalog cohorts. HGMD is unavailable and excluded. All 24 CADD GRCh37-v1.4 `inclAnno` batches were returned and audited.

## Files

- `analyze_cadd.py`, `audit_cadd_batches.py`, `collect_cadd_annotations.py`, and the two plotting scripts are the maintained code.
- `inputs/eqtl_gwas_variants_for_cadd_inclAnno.vcf.gz` is the submitted cohort VCF. It is a generated input, not a publisher source file.
- `outputs/` contains the 24 CADD returns, complete per-variant scores, audit/coverage tables, manifests, and final figures.
- The full method, metrics, and limitations are in [`FULL_ANALYSIS.md`](FULL_ANALYSIS.md).

## Remote copies

- **GitHub (code and documentation):** [ReproduceDeepsea/experiment/experiment9](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/experiment9)
- **Hugging Face (generated input and complete results):** [experiment/experiment9](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/experiment9)

The large `inputs/` and `outputs/` files are stored on Hugging Face; they are not duplicated in the GitHub repository. To restore them locally, use `hf download aer0vane/reproduce_deepsea --repo-type model --include 'experiment/experiment9/inputs/**' 'experiment/experiment9/outputs/**' --local-dir .` from the repository root.

## Reproduction boundary

The score pipeline reuses pretrained and locally retrained DeepSEA effect matrices from Figure 3; it does not run model inference or train a classifier. It calculates empirical E-values against each cohort's random-negative SNPs. See [`FULL_ANALYSIS.md`](FULL_ANALYSIS.md) for the formula, coverage, AUCs, and limitations.
