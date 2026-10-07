# Experiment 9 full analysis

Analysis used all 24 of 24 completed CADD GRCh37-v1.4 `inclAnno` batches. The scored files contain all eQTL and GWAS evaluation rows, using existing pretrained and ours DeepSEA effect matrices.

## Data audit and coverage

The CADD input VCF contains 1,855,581 unique SNVs. The returned files contain 1,855,578 unique variants after deduplication by CHROM/POS/REF/ALT. Every returned variant matches the submitted VCF; 3 submitted variants have no returned record. Overlapping batches have zero conflicting values for the four required fields. The Figure3 effect universe contains 1,855,704 variants, 123 of which were not included in the CADD VCF. Three returned rows have no usable value for any of the four conservation fields.

| Cohort/group | Variants | Any CADD feature | All four features |
|---|---:|---:|---:|
| eQTL positives | 78,612 | 78,610 (100.0%) | 78,440 (99.8%) |
| eQTL random negatives | 955,315 | 955,315 (100.0%) | 939,183 (98.3%) |
| eQTL 360bp negatives | 57,456 | 57,456 (100.0%) | 57,166 (99.5%) |
| GWAS positives | 12,296 | 12,295 (99.99%) | 12,271 (99.8%) |
| GWAS random negatives | 962,908 | 962,905 (100.0%) | 946,438 (98.3%) |
| GWAS 360bp negatives | 11,947 | 11,947 (100.0%) | 11,906 (99.7%) |

Full per-group coverage is in [`annotation_coverage.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/annotation_coverage.tsv); hashes, batch row counts, overlaps, and matching checks are in [`cadd_downloaded_batch_audit.json`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/cadd_downloaded_batch_audit.json).

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

The full AUC table also reports the 31kbp, 6.3kbp, 710bp, and 360bp negative groups in [`unsupervised_auc.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/unsupervised_auc.tsv). The two requested single-score panels are in [`suppfig6_pretrained.png`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/suppfig6_pretrained.png) and [`suppfig6_ours.png`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/suppfig6_ours.png); their x-axis labels use the cohort-specific mean distances measured in the reused evaluation sets, and n is the number of positive complete cases.

## Interpretation limits

These are unsupervised ranking results, not classifier results. AUCs against random negatives are modest (0.56–0.60); the combined score is highest among these three scores for both cohorts and models. Complete-case AUCs can be biased if missing annotation values are systematic. This analysis reuses Figure3 inference artifacts and did not rerun either checkpoint. HGMD is outside the planned scope of this experiment.

## Outputs

- [`analyze_cadd.py`](analyze_cadd.py): scoring and metric pipeline.
- [`eqtl_unsupervised_pretrained.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/experiment9/outputs/eqtl_unsupervised_pretrained.tsv.gz), [`eqtl_unsupervised_ours.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/experiment9/outputs/eqtl_unsupervised_ours.tsv.gz), [`gwas_unsupervised_pretrained.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/experiment9/outputs/gwas_unsupervised_pretrained.tsv.gz), [`gwas_unsupervised_ours.tsv.gz`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/experiment9/outputs/gwas_unsupervised_ours.tsv.gz): per-row features, E-values, and scores.
- [`unsupervised_auc.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/unsupervised_auc.tsv), [`annotation_coverage.tsv`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/annotation_coverage.tsv): metric and annotation coverage tables.
- [`analysis_manifest.json`](https://github.com/anaerovane/ReproduceDeepsea/blob/main/experiment/experiment9/outputs/analysis_manifest.json): source paths, completion snapshot, formulas, and limitations.
