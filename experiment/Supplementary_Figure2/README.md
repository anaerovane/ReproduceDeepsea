# Supplementary Figure 2 — DeepSEA vs gkm-SVM

This experiment compares DeepSEA with gkm-SVM across 690 transcription-factor predictors. The final evaluation has 689 evaluable predictors and one predictor with no held-out positives, for which AUC is undefined. The pretrained and local-checkpoint results are separate.

## Layout

- Root: training, retry, inference/plotting, and packaging scripts.
- `results/`: final AUC tables, metadata, audit notes, and the two final figures.
- Large trained gkm-SVM models and matched-test DeepSEA prediction caches are stored under this experiment on Hugging Face. The shared Princeton input bundle and local model weights are restored using `../../DOWNLOAD.md` and the project model resolver.

## Rebuild

Download the official gkm-SVM source from [Beer Lab](https://www.beerlab.org/gkmsvm/), compile `gkmsvm_kernel`, `gkmsvm_train`, and `gkmsvm_classify` into this folder, restore the Princeton data as described in `../../DOWNLOAD.md`, then run:

```bash
python reproduce_suppfig2_all_tf.py
python retry_incomplete_suppfig2.py
python plot_suppfig2_pretrained_ours.py
```

`package_suppfig2_artifacts.py` updates the matching Hugging Face experiment folder when new gkm-SVM model artifacts have been generated.

## Results

- `results/gkm_auc_results_completed.csv`: consolidated 690-predictor table.
- `results/pretrained_ours_same_test_auc.csv`: DeepSEA/gkm-SVM comparison on the same held-out examples.
- `results/gkm_vs_deepsea_pretrained.{png,pdf}` and `results/gkm_vs_deepsea_ours.{png,pdf}`: final figures.
- `results/COMPLETION_AUDIT.md` and `results/LOCAL_MODEL_FIGURES.md`: evaluation details and scope.
