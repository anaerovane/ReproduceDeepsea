# Supplementary Figure 2 — DeepSEA vs gkm-SVM

Reproduction of the transcription-factor binding comparison: the left panel compares per-feature DeepSEA and gkm-SVM AUCs; the right panel shows AUC distributions for DeepSEA (1000bp) and gkm-SVM (1000bp/300bp). Figures are split into **Pretrained** and **Ours** DeepSEA checkpoints.

This is a local reproduction on the documented DeepSEA held-out data and local model checkpoints. The values are computed from this project’s runs; they are not copied from the paper’s figure.

## Links

- [Hugging Face model files, trained gkm-SVM models, and result artifacts](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure2)
- [GitHub project](https://github.com/anaerovane/ReproduceDeepsea)
- [Original DeepSEA paper](https://www.nature.com/articles/nmeth.3547)
- [Official DeepSEA data/download page](https://deepsea.princeton.edu/help/)

## Files

- `reproduce_suppfig2_all_tf.py`: train/evaluate both gkm-SVM context lengths for all TF profiles.
- `retry_incomplete_suppfig2.py`: derive retry targets from the full-run results and retry them.
- `plot_suppfig2_pretrained_ours.py`: infer each DeepSEA checkpoint on matching held-out rows and draw separate figures.
- `package_suppfig2_artifacts.py`: merge per-feature working outputs into the single archive after a run.
- `gkm_vs_deepsea_pretrained.png/.pdf` and `gkm_vs_deepsea_ours.png/.pdf`: final plots.
- `gkm_auc_results_completed.csv` and `pretrained_ours_same_test_auc.csv`: AUC tables.
- `test_indices.npz`: held-out indices for the matched comparison.
- `vendor_gkmsvm_source.tar.zst` and `LICENSE.gkmsvm`: gkm-SVM source and license.
- The trained per-feature model/support-vector/score files are linked from the Hugging Face folder above; they are too large for this GitHub repository.

## Reproduce

Set `DEEPSEA_FUXIAN_ROOT` to a directory containing the official DeepSEA training bundle under `deepsea_train/`, the local checkpoints under `models/` and `training_checkpoints/`, and `predictor_names.txt`. The checkpoints are already available in the linked Hugging Face repository. Install Python dependencies (`numpy`, `scipy`, `torch`, `scikit-learn`, `matplotlib`) and run from this folder:

```bash
python reproduce_suppfig2_all_tf.py
python retry_incomplete_suppfig2.py
python plot_suppfig2_pretrained_ours.py
python package_suppfig2_artifacts.py
```

The included gkm-SVM executables are Linux x86-64. For another platform, build from the included source archive or use the upstream source linked in `LOCAL_MODEL_FIGURES.md`.

## Evaluation notes

The final table has 690 TF features: 689 have defined AUCs and one has no positive examples in the forward test split, so ROC AUC is undefined. Retry features use the maximum available balanced clean training sample; low-positive rows are identified in the CSV. The exact inputs, checkpoint hashes, sampling details, and limitations are in `LOCAL_MODEL_FIGURES.md` and `COMPLETION_AUDIT.md`.
