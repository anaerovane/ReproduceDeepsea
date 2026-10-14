# Local pretrained and ours TF comparison figures

The two final images are generated from real local DeepSEA checkpoint inference and gkm-SVM scores on the same per-feature held-out test indices:

- `gkm_vs_deepsea_pretrained.png` / `.pdf`: pretrained checkpoint retrieved by `model_assets.py` from Hugging Face or Zenodo (MD5 verified).
- `gkm_vs_deepsea_ours.png` / `.pdf`: ours checkpoint `best_model_FINAL_EPOCH53.pth`, resolved from Hugging Face when no local `training/checkpoints/` copy exists.

At the latest generation, 689 TF features had evaluable gkm-SVM AUC rows. Both DeepSEA checkpoints were run on the same saved forward-test indices (227,418 unique sequences), with strand predictions averaged. This includes the 16 retried features and their retry-specific held-out indices. Each plotted feature's DeepSEA AUC uses that feature's gkm test indices and labels. The gkm curves use their 1kb and center-300bp AUCs. Per-feature indices are in the Hugging Face `experiment/Supplementary_Figure2/results/test_indices.npz` cache; paired AUCs and test counts are in `pretrained_ours_same_test_auc.csv`; model hashes and run settings are in `pretrained_ours_same_test_metadata.json`.

The per-feature gkm model files, score files, and original individual index files are preserved in the Hugging Face archive. To unpack them for detailed inspection, download the archive from the linked Hub folder. After rerunning training, authenticate with `hf auth login` and run `python ../package_suppfig2_artifacts.py` from the experiment root; it merges new outputs with the Hub archive, uploads and verifies the result, and removes the local archive and expanded working folders.

## Data and model sources

- **Reference figure and method:** Zhou & Troyanskaya, *Nature Methods* (2015), [paper](https://www.nature.com/articles/nmeth.3547); [Supplementary Information PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM644_ESM.pdf), Supplementary Figure 2 (PDF page 3).
- **Sequence and TF labels:** the authors' [official DeepSEA site and download page](https://deepsea.princeton.edu/help/) provides the [training-data bundle](https://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz), including `train.mat` and `test.mat`. The bundle documents ENCODE/Roadmap chromatin labels over hg19 sequence. The locally used `train_data.bin`, `train_labels.bin`, and `test.mat` are from this DeepSEA training-data package.
- **Pretrained DeepSEA weights:** [Hugging Face model file](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/models/deepsea_predict.pth) and upstream [Zenodo record 1466993](https://zenodo.org/records/1466993), also indexed as [Kipoi DeepSEA/predict](https://kipoi.org/models/DeepSEA/predict/). Expected MD5 is `89e640bf6bdbe1ff165f484d9796efc7`.
- **Ours weights:** trained locally from the DeepSEA training bundle and published at `aer0vane/reproduce_deepsea:training_checkpoints/best_model_FINAL_EPOCH53.pth`. Its SHA256 is `5979e5b48ef892b76f23aa3af4149787e8f7a0ef5f5bfef5e676997a6c1731bb`.
- **gkm-SVM software:** the local run invokes the classic `gkmsvm_kernel`, `gkmsvm_train`, and `gkmsvm_classify` executables. The upstream software and source download are published by the [Beer Lab gkm-SVM page](https://www.beerlab.org/gkmsvm/); the algorithm is described in the [gkmSVM paper](https://academic.oup.com/bioinformatics/article/32/14/2205/1743153). The experiment locally trains the models; gkm AUCs are not imported from an external results table.

## How the plotted values were made

For each TF feature, up to 2,000 clean positive and negative sequences are sampled from the local DeepSEA training split. Separate gkm-SVM models are trained on 1,000bp and centered 300bp sequences. AUCs are measured on balanced examples from the forward half of the local DeepSEA test split, with up to 5,000 positive and the same number of negative examples. The 16 retry features use saved balanced training samples and test examples; their counts are recorded in `gkm_auc_results_completed.csv`. DeepSEA scores are freshly inferred from the local pretrained and locally trained checkpoints on the matched per-feature test sequences. The final plot data are local computations based on per-profile test indices.

Inputs use the verified `AGCT` channel order from the local test data for both checkpoints. The pretrained checkpoint uses its released `Threshold(0, 1e-6)` activation behavior; ours uses its saved epoch-53 architecture.
