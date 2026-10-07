# Experiment 2 — Context-length ablation

Compares pretrained and epoch-53 DeepSEA across centered 200, 500, and 1,000 bp input windows. The model weights are resolved through the shared `model_assets.py`; the test MAT file is restored from the official Princeton DeepSEA training bundle using [`../../DOWNLOAD.md`](../../DOWNLOAD.md).

## Files

- `experiment2_context_ablation.py`: computes per-predictor AUCs from the test set and writes `experiment2_auc_data.npz`.
- `experiment2_plot_overall.py`: draws one boxplot per model from that NPZ.
- `experiment2_auc_data.npz`: compact per-predictor results used by the plot script.
- `experiment2_overall_auc_pretrained.png` and `experiment2_overall_auc_ours.png`: final separate panels.

Code, the small result NPZ, and final figures are kept together in this experiment folder and published in the matching [GitHub experiment directory](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Experiment2_Context_Length_Ablation). These compact files stay on GitHub; the matching [Hugging Face experiment folder](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Experiment2_Context_Length_Ablation) is a navigation pointer only.
