# Supplementary Figure 4 — histone QTL effect prediction

This folder contains the current pretrained and locally retrained DeepSEA diagnostic figures for the histone-QTL task. The result is **not a full reproduction of the published curve**: the authors' original per-variant McVicker QTL/effect table was not recovered, so the cohort is reconstructed from GEO GSE47991 processed read counts.

## Contents

- `results/`: final pretrained/ours figures, compact accuracy table, and provenance. Raw predictions and per-variant tables are kept on Hugging Face.
- `inputs/`: the small GEO download manifest. Large public GEO tracks and the UCSC chain file are downloaded on demand by `download_inputs.py`.
- Root scripts: the CHT fitting, cohort summarization, inference, and plotting entry points.

## Data and models

- GEO processed H3K4me3/H3K27ac read-count tracks: [GSE47991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47991), 20 sample files enumerated with their direct NCBI FTP URLs and byte counts in `inputs/GEO_GSE47991_download_manifest.csv`.
- Coordinate conversion: UCSC [hg18-to-hg19 chain directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/); `download_inputs.py` retrieves `hg18ToHg19.over.chain.gz`.
- hg19 reference FASTA: use the shared project download instructions in `../../DOWNLOAD.md`.
- CHT implementation: [WASP](https://github.com/bmvdgeijn/WASP). Set `WASP_CHT_DIR` to a checkout containing its CHT scripts.
- The pretrained and local checkpoints are resolved by the shared `model_assets.py`; the recorded SHA-256 values are in `results/provenance.json`.
- Large raw prediction arrays and per-variant prediction tables are published under the matching [Hugging Face experiment folder](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure4/results); code, final plots, compact metrics, and this guide are on [GitHub](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Supplementary_Figure4). Restore the large files from the fuxian3 project root with `hf download aer0vane/reproduce_deepsea --include "experiment/Supplementary_Figure4/results/**" --local-dir .`.

## Rebuild

From this folder:

```bash
python download_inputs.py
export WASP_CHT_DIR=/path/to/WASP/CHT
bash run_mcvicker_cht.sh
```

The CHT rerun is computationally intensive. `run_mcvicker_cht.sh` regenerates the CHT tables, variant predictions, and figures under `results/`.

## Scope

The retained result uses the corrected hg18→hg19 mapping, reference-allele reconciliation, official DeepSEA 1000 bp VCF window convention, and forward/reverse-complement averaging for the pretrained variant-effects model. It contains 2,795 harmonized reconstructed SNP/mark rows (1,179 H3K4me3; 1,616 H3K27ac). The predictions and plots are diagnostic until the original McVicker per-variant association/effect cohort and exact analysis settings are available.
