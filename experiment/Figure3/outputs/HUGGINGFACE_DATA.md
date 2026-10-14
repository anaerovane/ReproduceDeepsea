# Figure 3 result files on Hugging Face

Browse or download the published outputs:

- [Hugging Face Figure 3 outputs](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs) (published revision `a01385444a8e82b59a33edfe26e62c6e6d4c5016`)
- [GitHub Figure 3 code, plots, and documentation](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)

Restore the published files to the `fuxian3` project root:

```bash
hf download aer0vane/reproduce_deepsea --include "experiment/Figure3/outputs/**" --local-dir .
```

## Uploaded files

The 28 files total about 27.48 GB. The two core model outputs are `ours_effects.npy` and `pretrained_effects.npy` (13,643,135,936 bytes each), containing ref/alt DeepSEA inference results. `variants.tsv` preserves the row order used by both matrices. The folder also contains four out-of-fold prediction tables, extracted CADD/GWAVA/FunSeq2 score tables, final AUC tables, and run manifests.

Pretrained and locally trained model weights are available from the [Hugging Face model repository](https://huggingface.co/aer0vane/reproduce_deepsea).
