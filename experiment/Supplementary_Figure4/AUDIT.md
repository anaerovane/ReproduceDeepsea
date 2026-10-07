# Supplementary Figure 4 result audit

## Retained output

`results/` contains one compact result set; the raw model-prediction NPZ and per-variant prediction tables are stored in the matching Hugging Face `experiment/Supplementary_Figure4/results/` folder. The retained result was generated using the hg18→hg19 chain in the correct direction, direct/complemented allele reconciliation against the hg19 reference, the DeepSEA VCF window centered at index 499, and forward/reverse-complement averaged pretrained variant-effect predictions.

The reconstructed cohort contains 2,795 harmonized SNP/mark rows: 1,179 H3K4me3 and 1,616 H3K27ac. Exact input/model hashes and inference settings are recorded in `results/provenance.json`.

## Reproduction boundary

These files are diagnostic results, not a complete reproduction of the paper's curve. The local cohort was recomputed from GEO GSE47991 donor read counts with WASP CHT and per-mark Benjamini–Hochberg FDR < 0.1. The authors' original per-variant McVicker QTL/effect table and exact cohort definition were not recovered. This can change cohort membership and effect directions.

The previous exploratory plots, partial chromosome runs, and runs with incorrect coordinate/allele handling were removed from the working folder. Do not use their old AUCs or figures.
