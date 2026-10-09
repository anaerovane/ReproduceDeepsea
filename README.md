# ReproduceDeepsea

Code, figures, evaluation tables, and selected large artifacts for reproducing DeepSEA results from Zhou & Troyanskaya, *Nature Methods* (2015).

- Paper: [Predicting effects of noncoding variants with deep learning–based sequence model](https://www.nature.com/articles/nmeth.3547)
- Large model weights and prediction artifacts: [Hugging Face repository](https://huggingface.co/aer0vane/reproduce_deepsea)
- Shared source-data and weight instructions: [`DOWNLOAD.md`](DOWNLOAD.md)

## Reproduction status

The project contains separate `pretrained` (released DeepSEA weights) and `ours` (locally retrained epoch-53 checkpoint) results. Completed figures do not imply that every result from the paper has been reproduced. HGMD positive variants are unavailable, so Figure 3 and Supplementary Figures 6–7 cover GRASP eQTL and GWAS Catalog cohorts only; ClinVar is not substituted for HGMD. Supplementary Figure 4 is a diagnostic reconstruction, and the HGMD indel task has not been completed.

## Experiments

Each experiment has its own folder under [`experiment/`](experiment/). The table links to the experiment guide in GitHub and to its corresponding Hugging Face directory. HF stores large data where needed; some HF folders are navigation pointers and are identified below.

| Figure or task | What is reproduced | GitHub code, figures, and notes | Hugging Face artifacts |
|---|---|---|---|
| Figure 2 | 919 chromatin-profile test predictions and DGF reference/alternative allele analysis; pretrained and ours plots | [Figure 2 folder](experiment/Figure2/README.md); shared inference/plot entrypoints are [`inference.py`](inference.py) and [`plot.py`](plot.py) | [Prediction arrays](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure2/outputs) |
| Figure 3 | GRASP eQTL and GWAS variant prioritization; DeepSEA compared with CADD, GWAVA, and FunSeq2; no HGMD | [Figure 3 guide, code, figures, and compact results](experiment/Figure3/README.md) | [Variant predictions, scores, and audit outputs](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs) |
| Supplementary Figure 1 | 200, 500, and 1,000 bp input-window comparison for pretrained and ours | [Code, AUC cache, and figures](experiment/Supplementary_Figure1/README.md) | [Folder](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure1) is a README pointer; result files are in GitHub |
| Supplementary Figure 2 | DeepSEA versus gkm-SVM on 690 TF profiles; 689 have defined held-out AUC, one has no test positives | [Code, final plots, and completion audit](experiment/Supplementary_Figure2/README.md) | [gkm-SVM models, indices, and prediction caches](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure2) |
| Supplementary Figure 3 | *In silico* saturation scans at four reported loci; 3,000 substitutions per locus, separate pretrained and ours plots | [Method, code, audit, and figures](experiment/Supplementary_Figure3/README.md) | [Sequences and prediction/effect caches](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure3) |
| Supplementary Figure 4 | Histone-QTL diagnostic reconstructed from GEO GSE47991 | [Code, provenance, and diagnostic plots](experiment/Supplementary_Figure4/README.md) | [Raw predictions and per-variant outputs](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure4/results) |
| Supplementary Figure 6 | Unsupervised functional-significance scores for GRASP and GWAS, across random and distance-based negative groups | [Code, full analysis, AUC tables, and pretrained/ours plots](experiment/Supplementary_Figure6/README.md) | [VCF, CADD responses, and per-variant scores](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure6) |
| Supplementary Figure 7 | Chromatin-only, conservation-only, and combined feature comparisons on GRASP/GWAS | [Code, method, and pretrained/ours plots](experiment/Supplementary_Figure7/README.md) | [OOF tables and fold caches](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Supplementary_Figure7/results) |
| Experiment 7 rerun | Direction check for three known FOXA1, GATA1, and FOXA2 variant mechanisms; 3/3 directions agree | [Script, raw predictions, and result note](experiment/experiment7_rerun/README.md) | [README pointer](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/experiment7_rerun); small outputs are in GitHub |

There is no `Supplementary_Figure5/` result folder in this project. The HGMD noncoding-indel experiment is also not completed.

## Important result boundaries

### Figure 2 code and data locations

The Figure 2 child folder contains its experiment README and input manifest. The inference/plot helpers are at the repository root, and prediction arrays are stored in the Hugging Face `experiment/Figure2/outputs/` folder. Running model inference also requires the DeepSEA runners, weights, test MAT file, Supplementary Table 3, and hg19 reference genome described in [`DOWNLOAD.md`](DOWNLOAD.md). The original allele-count table is [Springer Nature Supplementary Table 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM647_ESM.csv).

From the repository root, with the full DeepSEA project available at `PROJECT_ROOT`:

```bash
python inference.py --project-root "$PROJECT_ROOT" --device cuda
python plot.py pretrained
python plot.py ours
```

To redraw using published prediction caches, skip inference and run the two `plot.py` commands.

### Figure 3 scope and sampling

The DeepSEA-only OOF AUC table uses the full random-negative cohorts: 955,315 eQTL negatives and 962,908 GWAS negatives. The pretrained random-negative AUCs are 0.685 (eQTL) and 0.644 (GWAS); the corresponding ours values are 0.679 and 0.641. Results are chromatin-effects-only and omit the original model's four conservation features.

The external-method comparison table has only 47,336 eQTL and 47,623 GWAS random-negative common cases, about 4.95% of the full random-negative cohorts. [`baseline_source_manifest.json`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Figure3/outputs/baseline_source_manifest.json) records a 5% random-negative fraction (seed 42), while [`metrics_provenance.json`](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/experiment/Figure3/outputs/metrics_provenance.json) says `random_negative_cap=0`. Those records conflict. The external-method comparison should therefore be described as using the common scored subset, not as a full-random-negative comparison, until that provenance discrepancy is reconciled. CADD's version is unverified. Cohort overlap and allele/ref checks are documented in the Figure 3 provenance files.

To rerun from the Figure 3 directory after restoring the required data and weights:

```bash
python -u code/figure3_rerun.py > outputs/figure3_run.log 2>&1
```

### Other limitations

- Supplementary Figure 4 is a diagnostic: the original McVicker per-variant QTL/effect cohort was not recovered. Its reconstructed GSE47991 results are not the original paper curve.
- Supplementary Figure 6 evaluates the available GRASP and GWAS cohorts; HGMD is outside scope. The combined pretrained random-negative AUCs are 0.595 (eQTL) and 0.601 (GWAS); ours gives 0.589 and 0.598.
- Supplementary Figure 7 also uses GRASP/GWAS only and compares chromatin, conservation, and combined features. It does not include HGMD.
- Experiment 7 rerun is a three-locus direction check, not a full cohort reproduction.

## Reproducing individual experiments

Each linked experiment README lists required files, sources, commands, outputs, and limitations. Most large inputs can be restored selectively instead of downloading the whole Hugging Face repository. For example, run from the repository root:

```bash
hf download aer0vane/reproduce_deepsea \
  --include 'experiment/Figure3/outputs/**' \
  --local-dir .
```

Replace `Figure3` with the requested experiment folder. Figure 2 needs `experiment/Figure2/outputs/**`; Supplementary Figure 2 and 7 have their own large model/result subfolders. Supplementary Figure 1 and Experiment 7 rerun have small GitHub results and do not need large HF downloads.

For model weights, shared data, and reference-genome sources, see [`DOWNLOAD.md`](DOWNLOAD.md). Training scripts and training records are maintained in the local workspace; the published epoch-53 checkpoint is available from [Hugging Face](https://huggingface.co/aer0vane/reproduce_deepsea).

## Task B: separate external-method assessment

The local workspace also contains a separate `taskB/` analysis comparing DeepSEA with FunSeq2 and Nucleotide Transformer. It is not among the published `experiment/` folders in this repository. FunSeq2 plots/results currently reside in the local workspace. The NT-Multispecies 2.5B run was paused after 4,000 of 1,855,581 variants on CPU and has no full-cohort AUC; it requires a CUDA-visible GPU to resume. The NT paper's reported values are not local reproduction results.
