# Supplementary Figure 3: *In silico* saturated mutagenesis analysis for identifying predictive sequence features

From [Predicting effects of noncoding variants with deep learning–based sequence model](https://www-nature-com.libproxy1.nus.edu.sg/articles/nmeth.3547).

- [Original Supplementary Figure 3 image](https://media-springernature-com.libproxy1.nus.edu.sg/full/springer-static/image/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_Article_BFnmeth3547_Fig6_ESM.jpg)
- [GitHub source, code, and figures](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure3)
- [Hugging Face files and model checkpoints](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure3)

## Figure description

Predictive sequence features are identified by computationally mutating each base and measuring the effect on binding probability. The log2 fold change in odds is calculated from probability as *P* / (1 − *P*); yellow indicates increased binding probability and blue indicates decreased binding probability. The original figure shows whole-sequence scans and a 200 bp center view for four example sequences around chr1:109817590 G>T, chr16:209709 T>C, chr10:23508363 A>G, and chr16:52599188 C>T. The paper notes motifs including TTGCTCAA for CEBPB, TGATAA for GATA1, GTAAATA for FOXA1, and GTACATA for FOXA2.

## This reproduction package

This package runs the same four 1,000 bp sequence examples with the pretrained and locally trained DeepSEA checkpoints, then writes separate figures for each model. Results are model predictions, not experimental measurements. See [AUDIT.md](AUDIT.md) for input checks, effect calculations, and reproduction limits.

The experiment parameters are collected in `experiment6_config.json`. The four one-hot inputs are in `sequences/`. The sequence preparation script can regenerate them from GRCh37/hg19 chromosome FASTAs placed in the project-root `reference_genome/` directory. Large prediction and effect caches are stored under `results/` on Hugging Face, keeping this source folder limited to code, inputs, and final plots.

## Run from the project root

Use the `deepsea/fuxian3/` directory in this workspace, or the repository root in a GitHub checkout.

Download the shared feature-name table and checkpoints from Hugging Face into the matching paths at the repository root:

```bash
hf download aer0vane/reproduce_deepsea --local-dir . --include predictor_names.txt --include models/deepsea_predict.pth
```

The epoch-53 ours checkpoint is resolved from Hugging Face on demand by `model_assets.py`; it does not need to be copied into the project directory.

Then run inference and generate both figures:

```bash
python experiment/Supplementary_Figure3/experiment6_saturation_full.py
```

To redraw figures from the local effect cache only:

```bash
python experiment/Supplementary_Figure3/experiment6_plot.py
```

The raw prediction cache can be regenerated from the checkpoints and inputs. The final effect cache can be downloaded from [Hugging Face results](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure3/results) by `experiment6_plot.py` when no local cache exists. Both caches are retained there for audit and redraws.
