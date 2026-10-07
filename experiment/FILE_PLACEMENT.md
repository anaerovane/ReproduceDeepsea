# Experiment artifact placement

Checked on 2026-10-07. This table is the source of truth for artifact ownership. A remote copy is called verified only after checking its file list and sizes against local files. Large local copies are removed only after that check.

## Placement rules

- **GitHub:** training/inference/plotting code, README and audit/provenance documents, final figures, and small result tables needed to understand the figure.
- **Hugging Face:** model checkpoints, large generated predictions/fold caches, and returned annotations that are expensive or impossible to regenerate exactly. Each experiment has its own named path.
- **Local only:** active work caches that are not yet safely archived, or temporary data needed to continue a not-yet-finished task. Publisher/source data should be restored from its official URL in `../DOWNLOAD.md` rather than kept as another project copy.

## Verified and pending locations

| Experiment | GitHub | Hugging Face | Local |
|---|---|---|---|
| Figure 2 / DGF | `experiment/Figure2`: code, README, manifest, final figures | `experiment/Figure2/outputs`: generated inference arrays and test labels. The publisher source CSV is linked from Springer Nature, not mirrored. | code/docs/figures only; large copies removed after remote verification |
| Context-length ablation (`Supple1`) | Not yet published | No project result files verified | Current small code, two figures, and NPZ remain local; upload/README cleanup pending |
| Supplementary Figure 2 / gkm-SVM | `experiment/Supplementary_Figure2`: scripts, docs, final pretrained/ours figures and compact AUC tables | `experiment/Supplementary_Figure2/feature_models_and_indices.tar.zst`, `test_indices.npz`, and run metadata | Working code and flat prediction/test-input caches remain; still being audited before pruning |
| Histone QTL (`experiment5_histone_qtl`) | Not yet published | No final artifact set verified | Large historical/source/result tree remains; only the corrected-coordinate diagnostic is valid, not a paper reproduction. Cleanup pending |
| Figure 3 / variant prioritization | `experiment/Figure3`: code, audit/provenance, final figures and guide | `experiment/Figure3/outputs`: large effects, variant mapping, OOF and baseline result artifacts | Code/docs/figures; large duplicate caches removed where verified |
| Supplementary Figure 3 / saturation | `experiment/Supplementary_Figure3`: code, audit, final pretrained/ours figures | `experiment/Supplementary_Figure3/results`: raw prediction and derived effect NPZs plus metadata | Code/docs/figures and the small sequence inputs |
| Known-variant mechanism check (`experiment7_rerun`) | Not yet published | No remote copy verified | Script, audit JSON, raw NPZ, and report remain local; upload/README cleanup pending |
| Figure 6 / unsupervised significance (`experiment9`) | `experiment/experiment9`: maintained scripts and full analysis docs | `experiment/experiment9/inputs` and `outputs`: submitted VCF, all 24 CADD returns, complete score tables, audits and figures | Code/docs only; verified generated/input files removed locally |
| Supplementary Figure 7 / feature ablation | `experiment/Supplementary_Figure7`: code, README, final pretrained/ours figures | `experiment/Supplementary_Figure7/results`: 112 fold/OOF/summary files | Code/docs/final figures only; verified result cache removed locally |

## Important boundaries

- Figure 3, Figure 6, and Supplementary Figure 7 omit HGMD because its licensed professional dataset was unavailable. ClinVar is not presented as HGMD.
- Figure 6's 24/24 CADD results and Supplementary Figure 7's 112 result files have been checked on Hugging Face by exact filename and size. Figure 7 images were checked on GitHub by exact filename and size.
- Experiment 3 is complete at the evaluation-table level: 689 predictors have AUCs; one predictor has no positive test examples and is unscorable. Its remaining local large-cache audit is still pending.
- Experiment 5 remains diagnostic only; do not label its corrected-coordinate plot as the paper result.
