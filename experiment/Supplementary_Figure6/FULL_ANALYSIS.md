# Supplementary Figure 6 full analysis

Analysis used all 24 of 24 completed CADD GRCh37-v1.4 `inclAnno` batches. The scored files contain all eQTL and GWAS evaluation rows, using existing pretrained and ours DeepSEA effect matrices.

## Method

For each of 919 chromatin tracks, the effect magnitude is `|P_ref - P_alt| * |logit(P_ref) - logit(P_alt)|`. Each track's empirical E-value is the proportion of that cohort's random-negative background with a strictly larger value, using add-one smoothing. The four CADD fields (`priPhCons`, `priPhyloP`, `GerpN`, `GerpS`) are separately converted to empirical upper-tail E-values against the same random-negative background. The chromatin and conservation scores are geometric means of their feature E-values; the combined functional-significance score is their product. Lower scores indicate stronger functional significance.

AUCs below use the negative of each score and compare positives with random-negative SNPs. Chromatin AUCs use all cohort rows; conservation and combined AUCs use rows complete for all four CADD fields.

## AUC against random-negative SNPs

| Cohort | Model | Score | Positives | Negatives | AUC |
|---|---|---|---:|---:|---:|
| eQTL | pretrained | Chromatin | 78,612 | 955,315 | 0.580 |
| eQTL | pretrained | Conservation | 78,440 | 939,183 | 0.561 |
| eQTL | pretrained | Combined | 78,440 | 939,183 | 0.595 |
| eQTL | ours | Chromatin | 78,612 | 955,315 | 0.570 |
| eQTL | ours | Conservation | 78,440 | 939,183 | 0.561 |
| eQTL | ours | Combined | 78,440 | 939,183 | 0.589 |
| GWAS | pretrained | Chromatin | 12,296 | 962,908 | 0.567 |
| GWAS | pretrained | Conservation | 12,271 | 946,438 | 0.579 |
| GWAS | pretrained | Combined | 12,271 | 946,438 | 0.601 |
| GWAS | ours | Chromatin | 12,296 | 962,908 | 0.563 |
| GWAS | ours | Conservation | 12,271 | 946,438 | 0.579 |
| GWAS | ours | Combined | 12,271 | 946,438 | 0.598 |

The full AUC table also reports the 31kbp, 6.3kbp, 710bp, and 360bp negative groups in [`unsupervised_auc.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/unsupervised_auc.tsv). The two requested single-score panels are in [`Supplementary_Figure6_pretrained.png`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/Supplementary_Figure6_pretrained.png) and [`Supplementary_Figure6_ours.png`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/Supplementary_Figure6_ours.png); their x-axis labels use the cohort-specific mean distances measured in the reused evaluation sets, and n is the number of positive complete cases.

## Outputs

- [`analyze_cadd.py`](analyze_cadd.py): scoring and metric pipeline.
- [`eqtl_unsupervised_pretrained.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Supplementary_Figure6/outputs/eqtl_unsupervised_pretrained.tsv.gz), [`eqtl_unsupervised_ours.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Supplementary_Figure6/outputs/eqtl_unsupervised_ours.tsv.gz), [`gwas_unsupervised_pretrained.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Supplementary_Figure6/outputs/gwas_unsupervised_pretrained.tsv.gz), [`gwas_unsupervised_ours.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Supplementary_Figure6/outputs/gwas_unsupervised_ours.tsv.gz): per-row features, E-values, and scores.
- [`unsupervised_auc.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/unsupervised_auc.tsv), [`annotation_coverage.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/annotation_coverage.tsv): metric and annotation coverage tables.
- [`analysis_manifest.json`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/Supplementary_Figure6/outputs/analysis_manifest.json): source paths, completion snapshot, and formulas.
