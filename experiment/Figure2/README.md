# DeepSEA Figure 2

This folder contains the inference orchestrator, plotting code, manifest, and the two generated plots. Prediction arrays are kept separately under `outputs/`; the published experimental table belongs under `inputs/`.

## Contents

- `inference.py`: runs the pretrained and ours inference entry points; fresh arrays go to `outputs/`.
- `plot.py`: generates one plot for `pretrained` or `ours`, using local inputs or downloading the matching Hugging Face files.
- `figure2_pretrained.png`, `figure2_ours.png`: separate generated plots.
- `input_manifest.json`: file sizes, checksums, and source paths for the inputs used for plotting.
- `inputs/41592_2015_BFnmeth3547_MOESM647_ESM.csv`: experimental allele count data from the paper's Supplementary Table 3.
- `outputs/test_labels.npy`: labels for the DeepSEA test set.
- `outputs/test_predictions_{pretrained,ours}.npy`: model scores on the test set.
- `outputs/dgf_{ref,alt}_preds_{pretrained,ours}.npy`: scores for reference and alternative DGF sequences.
- `../../predictor_names.txt`: shared labels for the 919 chromatin output features.

## Recreate the plots

From this folder, with the fuxian3 model and inference inputs available:

```bash
python inference.py --device cuda
python plot.py pretrained
python plot.py ours
```

`inference.py` calls `../../run_inference.py` and `../../run_inference_ours.py`. It requires `deepsea_train/test.mat`, the Supplementary Table 3 CSV under `newdata/`, and the reference genome under `reference_genome/`. Model weights are resolved by `model_assets.py`. Restore downloadable inputs using [`DOWNLOAD.md`](../../DOWNLOAD.md).

To redraw from the published prediction caches without rerunning inference:

```bash
python plot.py pretrained
python plot.py ours
```

The plotting script downloads the NPY caches from [`Hugging Face: Figure 2 outputs`](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure2/outputs) as needed. `predictor_names.txt` is a shared root-level metadata file. The allele-count table is Supplementary Table 3 from Zhou & Troyanskaya, *Nature Methods* (2015), available directly from [Springer Nature](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM647_ESM.csv).

The `.npy` arrays are generated inference outputs, not original paper downloads. The PNGs are our separate pretrained and ours reproductions, not copies of the paper's image. There are no `.npz` inputs for this figure.
