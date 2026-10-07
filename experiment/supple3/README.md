# Experiment 6: saturation mutagenesis

This folder contains the configuration, sequence preparation, inference and plotting code, audit notes, and final pretrained/ours figures for the saturation mutagenesis analysis.

## What this experiment does

For four configured 1,000 bp GRCh37/hg19 reference windows, it evaluates all 3,000 single-base substitutions per window with each of two existing DeepSEA checkpoints. It saves model probabilities and log2 odds effects, then plots each model separately. **This experiment does not train either checkpoint.** The checkpoints are shared project artifacts produced or obtained outside this experiment; their training and provenance belong in the project-level model documentation.

The effect is `log2(mutant odds / reference odds)`. These are in-silico model predictions, not experimental measurements. See [AUDIT.md](AUDIT.md) for the input checks and limitations.

## Inputs

The single source of experiment parameters is [`experiment6_config.json`](experiment6_config.json). It declares the four variants, reference and alternate alleles, TF/cell labels, output-feature index and expected name, sequence channel order, checkpoint paths, and plotting window. The plotting code validates every configured feature index against `predictor_names.txt` and stops on a mismatch; it does not infer or silently substitute a feature.

The four 1 kb one-hot sequence inputs are in the single-level `sequences/` folder. Regenerate them from `deepsea/fuxian3/reference_genome/{chrom}.fa` with:

```bash
python experiment/supple3/prepare_sequences.py
```

The preparation script uses 1-based GRCh37/hg19 coordinates, checks the reference allele, rejects ambiguous bases, and writes the one-hot arrays. It generates the same four arrays as the inputs used for the final inference.

Download the shared feature-name table and checkpoints from Hugging Face into these paths at the repository root:

```bash
hf download aer0vane/reproduce_deepsea --local-dir . --include predictor_names.txt --include models/deepsea_predict.pth --include training_checkpoints/best_model_FINAL_EPOCH53.pth
``` Checkpoint architectures and channel/reverse-complement handling match the corresponding model definitions in the inference script and are loaded strictly.

The configured samples are chr1:109817590 G>T, chr16:209709 T>C, chr10:23508363 A>G, and chr16:52599188 C>T. The alternate allele documents the variant; saturation mutagenesis evaluates all three alternate bases at every sequence position.

## Run

From the repository root, after the project inputs and checkpoints are present:

```bash
python experiment/supple3/experiment6_saturation_full.py
```

This runs both checkpoints, rewrites the two local `.npz` caches, and generates the figures. To redraw only from the existing effect cache:

```bash
python experiment/supple3/experiment6_plot.py
```

The environment needs Python with NumPy, PyTorch, and Matplotlib. The full inference script uses CUDA when available and otherwise runs on CPU. The verified run used the project's `deepsea` conda environment.

## Tracked outputs

- `experiment6_saturation_pretrained.png` and `.pdf`
- `experiment6_saturation_ours.png` and `.pdf`
- `experiment6_plot_metadata.json`: selected output-feature indices and plotting limits
- `AUDIT.md`: checks, interpretation, and reproduction boundary

The `.npz` prediction/effect caches are retained locally for provenance and redraws, but excluded from ordinary Git commits because together they are about 163 MiB. Keep them in a separate artifact store if they need to be shared; do not mistake them for model weights.

## Scope

The model architecture and checkpoint locations are fixed to the two named project models. Experiment values (variant coordinates, expected bases, output-feature indices, checkpoint paths, sequence channel order, and zoom interval) live in the JSON configuration rather than being duplicated between scripts. Calculated effects, color limits, and figure values come from the loaded checkpoints and inference caches; none are manually entered as plotted results.
