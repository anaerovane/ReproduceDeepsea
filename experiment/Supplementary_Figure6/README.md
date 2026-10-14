# Supplementary Figure 6 — Unsupervised functional significance

The completed analysis uses GRASP eQTL and GWAS Catalog cohorts with 24 CADD GRCh37-v1.4 `inclAnno` batches.

## Files

- `analyze_cadd.py`, `audit_cadd_batches.py`, `collect_cadd_annotations.py`, and the two plotting scripts are the maintained code.
- `inputs/eqtl_gwas_variants_for_cadd_inclAnno.vcf.gz`: cohort VCF used for CADD annotation.
- `outputs/` on Hugging Face contains the 24 CADD returns and complete per-variant score files. Compact audit tables, the AUC summary, and final Supplementary Figure 6 panels are on GitHub under `outputs/`.
- The method and metrics are in [`FULL_ANALYSIS.md`](FULL_ANALYSIS.md).

## Remote copies

- **GitHub (code and documentation):** [ReproduceDeepsea/experiment/Supplementary_Figure6](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure6)
- **GitHub (compact results and final panels):** [experiment/Supplementary_Figure6/outputs](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure6/outputs)
- **Hugging Face (generated input and full result files):** [experiment/Supplementary_Figure6](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure6)

The large `inputs/` and `outputs/` files are stored on Hugging Face. To restore them locally, use `hf download aer0vane/reproduce_deepsea --repo-type model --include 'experiment/Supplementary_Figure6/**' --local-dir .` from the repository root.

The score pipeline calculates empirical E-values from pretrained and locally retrained DeepSEA effect matrices and cohort random-negative SNPs. See [`FULL_ANALYSIS.md`](FULL_ANALYSIS.md) for formulas, coverage, and AUCs.
