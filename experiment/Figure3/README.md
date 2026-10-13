# Figure 3: Functional Variant Prioritization

This Figure 3 analysis compares pretrained and locally trained DeepSEA models with CADD PHRED, GWAVA (unmatched, TSS, and region), and FunSeq2 v2.1.0 on GRASP eQTL and GWAS Catalog cohorts.

For each variant, the pipeline constructs 1,000 bp hg19 reference and alternate sequences, validates the reference allele, and evaluates both DNA strands. DeepSEA scores summarize reference-to-alternate changes across 919 predicted chromatin profiles. Evaluation preserves the cohort fold assignments and calculates out-of-fold ROC AUC for each control-variant distance group. Sample counts, mean distances, and scores are derived from the processed variant records.

## Figures

- Pretrained DeepSEA: `figures/figure3_pretrained_separate.png`
- Locally trained DeepSEA: `figures/figure3_ours_separate.png`

## Code and results

- Reproduction code and figures: [GitHub Figure 3](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)
- Inference outputs and evaluation tables: [Hugging Face Figure 3 outputs](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs)
