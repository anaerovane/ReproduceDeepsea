# ReproduceDeepsea

Reproduction code, figures, evaluation tables, and model artifacts for DeepSEA, based on Zhou & Troyanskaya, *Nature Methods* (2015).

- Paper: [Predicting effects of noncoding variants with deep learning–based sequence model](https://www.nature.com/articles/nmeth.3547)
- Model weights and large prediction artifacts: [Hugging Face repository](https://huggingface.co/aer0vane/reproduce_deepsea)
- Shared data and weight download instructions: [`DOWNLOAD.md`](DOWNLOAD.md)

## Experiments

Each entry links to its code, figures, and experiment guide on GitHub, together with its matching Hugging Face artifact directory.

| Experiment | Contents | GitHub | Hugging Face |
|---|---|---|---|
| Figure 2 | Chromatin-profile test predictions and DGF reference/alternative allele analysis; pretrained and ours | [Code, plots, and guide](experiment/Figure2/README.md) | [Prediction arrays](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure2/outputs) |
| Figure 3 | GRASP eQTL and GWAS variant prioritization with DeepSEA, CADD, GWAVA, and FunSeq2 | [Code, plots, and results](experiment/Figure3/README.md) | [Variant scores and result artifacts](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs) |
| Supplementary Figure 1 | Input-window comparison at 200, 500, and 1,000 bp for pretrained and ours | [Code, plots, and results](experiment/Supplementary_Figure1/README.md) | [Experiment directory](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure1) |
| Supplementary Figure 2 | DeepSEA and gkm-SVM transcription-factor binding comparison | [Code, plots, and results](experiment/Supplementary_Figure2/README.md) | [Models and prediction artifacts](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure2) |
| Supplementary Figure 3 | *In silico* saturation scans at four reported loci for pretrained and ours | [Code, plots, and results](experiment/Supplementary_Figure3/README.md) | [Sequence and prediction artifacts](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure3) |
| Supplementary Figure 4 | Histone-QTL analysis using GSE47991 | [Code, plots, and results](experiment/Supplementary_Figure4/README.md) | [Prediction and variant results](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure4/results) |
| Supplementary Figure 6 | Functional-significance scores for GRASP eQTL and GWAS variants | [Code, plots, and results](experiment/Supplementary_Figure6/README.md) | [VCF, CADD scores, and variant results](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure6) |
| Supplementary Figure 7 | Chromatin, conservation, and combined-feature comparisons for GRASP and GWAS | [Code, plots, and results](experiment/Supplementary_Figure7/README.md) | [Fold and OOF results](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure7/results) |
| Experiment 7 rerun | Variant-effect direction analysis at FOXA1, GATA1, and FOXA2 loci | [Code, plots, and results](experiment/experiment7_rerun/README.md) | [Experiment directory](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/experiment7_rerun) |

## Figure 2 inference and plotting

The Figure 2 inference and plotting entry points are at the repository root. With the DeepSEA project files restored as described in [`DOWNLOAD.md`](DOWNLOAD.md), run:

```bash
python inference.py --project-root "$PROJECT_ROOT" --device cuda
python plot.py pretrained
python plot.py ours
```

To redraw the figures from prediction arrays, run the two `plot.py` commands. Inference arrays are published in the Hugging Face Figure 2 folder.

## Figure 3

From `experiment/Figure3/`, run the variant-prioritization analysis with:

```bash
python -u code/figure3_rerun.py
```

See [`experiment/Figure3/README.md`](experiment/Figure3/README.md) for inputs, data sources, and outputs. Large experiment files can be selectively downloaded from Hugging Face; for example:

```bash
hf download aer0vane/reproduce_deepsea \
  --include 'experiment/Figure3/outputs/**' \
  --local-dir .
```

Replace `Figure3` with the experiment directory to retrieve its artifacts. Shared model weights, source data, and reference genome instructions are in [`DOWNLOAD.md`](DOWNLOAD.md).
