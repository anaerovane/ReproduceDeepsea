# Supplementary Figure 1: Context-length ablation

Compares pretrained and epoch-53 DeepSEA across centered 200, 500, and 1,000 bp input windows. The model weights are resolved through the shared `model_assets.py`; the test MAT file is restored from the official Princeton DeepSEA training bundle using [`../../DOWNLOAD.md`](../../DOWNLOAD.md).

## Files

- `suppfig1_context_ablation.py`: computes per-predictor AUCs from the test set and writes `supplementary_figure1_auc_data.npz`.
- `plot_suppfig1_overall.py`: draws one boxplot per model from that NPZ.
- `supplementary_figure1_auc_data.npz`: compact per-predictor results used by the plot script.
- `supplementary_figure1_pretrained.png` and `supplementary_figure1_ours.png`: final separate panels.
- `supplementary_figure1_context_ablation.png`: feature-class comparison produced by the AUC computation script.

The code, result NPZ, and final figures are kept together in this directory and published in the matching [GitHub folder](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure1). The [Hugging Face folder](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure1) links back to GitHub.

Run from the `deepsea/` project root:

```bash
python experiment/Supplementary_Figure1/suppfig1_context_ablation.py
python experiment/Supplementary_Figure1/plot_suppfig1_overall.py
```
