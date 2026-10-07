# Local pretrained and ours TF comparison figures

The two final images are generated from real local DeepSEA checkpoint inference and gkm-SVM scores on the same per-feature held-out test indices:

- `gkm_vs_deepsea_pretrained.png` / `.pdf`: pretrained checkpoint `deepsea/fuxian3/models/deepsea_predict.pth`.
- `gkm_vs_deepsea_ours.png` / `.pdf`: ours checkpoint `deepsea/fuxian3/training_checkpoints/best_model_FINAL_EPOCH53.pth`.

At the latest generation, 689 TF features had evaluable gkm-SVM AUC rows. Both DeepSEA checkpoints were run on the union of the same saved forward-test indices (227,418 unique sequences), with strand predictions averaged. This includes the 16 retried features and their retry-specific held-out indices. The one remaining feature has no positive examples in the forward test split and therefore has no ROC AUC. Each plotted feature's DeepSEA AUC uses exactly that feature's gkm test indices and labels. The gkm curves use their 1kb and center-300bp AUCs. Per-feature indices are in `test_indices.npz`; per-feature paired AUCs and test counts are in `pretrained_ours_same_test_auc.csv`; model hashes and run settings are in `pretrained_ours_same_test_metadata.json`.

The per-feature gkm model files, score files, and original individual index files are preserved in the Hugging Face archive. To unpack them for detailed inspection, download the archive from the linked Hub folder. After rerunning training, authenticate with `hf auth login` and run `python ../package_suppfig2_artifacts.py` from the experiment root; it merges new outputs with the Hub archive, uploads and verifies the result, and removes the local archive and expanded working folders.

## Data and model sources

- **Reference figure and method:** Zhou & Troyanskaya, *Nature Methods* (2015), [paper](https://www.nature.com/articles/nmeth.3547); [Supplementary Information PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM644_ESM.pdf), Supplementary Figure 2 (PDF page 3).
- **Sequence and TF labels:** the authors' [official DeepSEA site and download page](https://deepsea.princeton.edu/help/) provides the [training-data bundle](https://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz), including `train.mat` and `test.mat`. The bundle documents ENCODE/Roadmap chromatin labels over hg19 sequence. The locally used `train_data.bin`, `train_labels.bin`, and `test.mat` are from this DeepSEA training-data package.
- **Pretrained DeepSEA weights:** [Zenodo record 1466993](https://zenodo.org/records/1466993), also indexed as [Kipoi DeepSEA/predict](https://kipoi.org/models/DeepSEA/predict/). The local `models/deepsea_predict.pth` has MD5 `89e640bf6bdbe1ff165f484d9796efc7`, matching the Zenodo file record.
- **Ours weights:** `training_checkpoints/best_model_FINAL_EPOCH53.pth` was trained locally from the DeepSEA training bundle; it is not a downloaded public checkpoint. Its SHA256 is `5979e5b48ef892b76f23aa3af4149787e8f7a0ef5f5bfef5e676997a6c1731bb`.
- **gkm-SVM software:** the local run invokes the classic `gkmsvm_kernel`, `gkmsvm_train`, and `gkmsvm_classify` executables. The upstream software and source download are published by the [Beer Lab gkm-SVM page](https://www.beerlab.org/gkmsvm/); the algorithm is described in the [gkmSVM paper](https://academic.oup.com/bioinformatics/article/32/14/2205/1743153). The experiment locally trains the models; gkm AUCs are not imported from an external results table.

## How the plotted values were made

For each TF feature, positive and negative sequences are sampled from the local DeepSEA training split (up to 2,000 clean sequences per class), then separate gkm-SVM models are trained on 1,000bp and centered 300bp sequences. AUCs are measured on a balanced subset of the forward half of the local DeepSEA test split (up to 5,000 positives and the same number of negatives). The 16 retry features use the available balanced clean training samples and all available test positives; their exact sample counts are recorded in `gkm_auc_results_completed.csv`. The plotted DeepSEA AUCs are freshly inferred from the local pretrained or locally trained checkpoint on those same per-feature test sequences. The 0-positive feature is omitted because ROC AUC is undefined.

The final plot data are local computations, not values copied from Supplementary Table 2. The older `gkm_vs_deepsea_approx.*` draft used published per-profile AUCs and is superseded.

The original test cohort is stratified per TF: up to 5,000 positive examples and an equal number of negative examples from the forward half of `test.mat`. Retry cohorts use all available positives and the same number of negatives. Consequently, different TFs can have different test-set sizes; a few retried TFs have fewer than 10 positives, so their AUCs are statistically unstable. This reproduces the requested comparison layout on the available local held-out data; it is not a claim that all paper-specific split and sampling details were recovered.

Inputs use the verified `AGCT` channel order from the local test data for both checkpoints. The pretrained checkpoint uses its released `Threshold(0, 1e-6)` activation behavior; ours uses its saved epoch-53 architecture. The previous `gkm_vs_deepsea_approx.*` draft used published Supplementary Table 2 values and is superseded by these two local-checkpoint figures.
