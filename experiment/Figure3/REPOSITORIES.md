# Figure 3 files

- GitHub: [code, plots, method notes, and compact result summaries](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)
- Hugging Face: [core inference outputs and evaluation artifacts](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs)

Hugging Face revision [`a01385444a8e82b59a33edfe26e62c6e6d4c5016`](https://huggingface.co/aer0vane/reproduce_deepsea/tree/a01385444a8e82b59a33edfe26e62c6e6d4c5016/experiment/Figure3/outputs) contains 28 selected files (about 27.48 GB): two model effect matrices, their variant-row mapping, OOF predictions, baseline-score tables, AUC summaries, and audit metadata. It excludes reconstructable sequence caches and fold-level classifier caches.

Restore the published files to the `fuxian3` project root:

```bash
hf download aer0vane/reproduce_deepsea --include "experiment/Figure3/outputs/**" --local-dir .
```

Model weights remain in the existing Hugging Face model paths; GitHub contains no model checkpoints or large inference arrays.
