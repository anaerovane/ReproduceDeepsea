# Experiment 3 — Supplementary Figure 2

Reproduction code for the DeepSEA vs gkm-SVM transcription-factor comparison. The local pretrained and locally trained checkpoints have separate figures.

- [GitHub code, figures, and documentation](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure2)
- [Hugging Face trained gkm-SVM models and result files](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure2)

## Directory layout

- This directory contains the training, retry, inference/plotting, and packaging scripts, gkm-SVM executables, and this guide.
- `suppfig2_tf_fulltrain/` is the only child directory. Its images, tables, model inference outputs, and notes are directly inside it.
- Detailed per-feature models and training indices are on Hugging Face at `experiment/Supplementary_Figure2/feature_models_and_indices.tar.zst`. DeepSEA matched-test prediction caches and index metadata are under `experiment/Supplementary_Figure2/results/`; the plotting script fetches those caches on demand. The large remote copies are removed from the local experiment folder after checksum verification.
- Upstream gkm-SVM source and superseded smoke/legacy run materials are retained in separate archives in that same child directory.
- The archives preserve files and their original paths; they can be expanded when per-feature inspection or compilation is needed.

## Active scripts

- `reproduce_suppfig2_all_tf.py` trains and evaluates the 1000bp and center-300bp gkm-SVM models.
- `retry_incomplete_suppfig2.py` retries features that lack a valid full-run AUC.
- `plot_suppfig2_pretrained_ours.py` evaluates the pretrained and local checkpoints on matching held-out examples and writes the two figures.

Set `DEEPSEA_FUXIAN_ROOT` to the `fuxian3` project root. Restore missing data using [`DOWNLOAD.md`](../../DOWNLOAD.md); model weights are fetched into the user cache by `model_assets.py`. Pretrained weights are downloaded on demand by `model_assets.py`; the epoch-53 checkpoint is resolved from Hugging Face if it is not present in `training/checkpoints/`. Run from this directory:

```bash
python prepare_gkm_training_cache.py  # optional: the runner creates caches automatically when absent
python reproduce_suppfig2_all_tf.py
python retry_incomplete_suppfig2.py
python plot_suppfig2_pretrained_ours.py
python package_suppfig2_artifacts.py
```

The training scripts create per-feature working folders while they run. Log in with `hf auth login` first; the final command temporarily downloads the existing Hub archive, merges new outputs, uploads and verifies the updated archive, then removes expanded folders and all temporary large copies. Set `DEEPSEA_HF_REPO` to override the default Hub repository.

## Final outputs

- `suppfig2_tf_fulltrain/gkm_vs_deepsea_pretrained.png` and `.pdf`
- `suppfig2_tf_fulltrain/gkm_vs_deepsea_ours.png` and `.pdf`
- `suppfig2_tf_fulltrain/gkm_auc_results_completed.csv`
- `suppfig2_tf_fulltrain/pretrained_ours_same_test_auc.csv`

The published-table proxy chart was removed; it was superseded by the two local-checkpoint comparisons. See `suppfig2_tf_fulltrain/LOCAL_MODEL_FIGURES.md` for methods, sources, cohort details, and limitations.
