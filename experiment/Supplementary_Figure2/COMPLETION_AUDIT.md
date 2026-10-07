# Supplementary Figure 2 gkm-SVM completion audit

- The original run produced AUCs for 673/690 TF features; 11 training-sample errors and 6 test-sample skips were retried.
- Sixteen retry features produced 1000bp and center-300bp gkm-SVM AUCs using the maximum available balanced clean training sample and every available forward-test positive matched to the same number of negatives.
- Feature 598 (`K562_BRF1_None`) has 0 positive labels in the forward test split. ROC AUC is undefined, so it remains explicitly unscorable.
- `gkm_auc_results_completed.csv` is the consolidated table: 689 evaluable features and 1 unscorable feature. Rows with fewer than 10 positives are statistically uncertain; the feature with 1 positive is a ranking against that single observation.
- The superseded initial-run CSV and transient run/smoke logs were removed to prevent confusion. The final consolidated CSV is the authoritative table; per-feature gkm models and source indices are in the Hugging Face archive.
- Per-feature held-out indices are available for direct use in `test_indices.npz`; this flat file prefers retry-specific indices where applicable.
- The final Pretrained and Ours figures use the consolidated 689-feature cohort and freshly inferred DeepSEA checkpoint scores on the matching held-out indices.
