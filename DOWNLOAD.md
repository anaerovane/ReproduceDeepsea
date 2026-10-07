# Data and model download guide

This is the single download guide for the DeepSEA reproduction workspace. Large source data and optional model files are not kept in this checkout. Run commands from the `fuxian3` project root unless a section says otherwise. The download targets below match the paths used by the existing scripts.

## 1. DeepSEA train, validation, and test data

The official Princeton bundle is a 3,821,872,019-byte gzip archive (about 3.6 GiB compressed). It contains the `.mat` data used for local training and inference. The current code expects `deepsea_train/train.mat`, `deepsea_train/valid.mat`, and `deepsea_train/test.mat`.

- Instructions: <https://deepsea.princeton.edu/help/>
- Direct bundle: <https://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz>

```bash
mkdir -p /tmp/deepsea_bundle_unpack deepsea_train
curl -fL --retry 3 \
  -o /tmp/deepsea_train_bundle.v0.9.tar.gz \
  https://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz
tar -xzf /tmp/deepsea_train_bundle.v0.9.tar.gz -C /tmp/deepsea_bundle_unpack
bundle_dir="$(dirname "$(find /tmp/deepsea_bundle_unpack -type f -name train.mat -print -quit)")"
cp "$bundle_dir"/train.mat "$bundle_dir"/valid.mat "$bundle_dir"/test.mat deepsea_train/
```

The MAT files encode 1 kb sequences and labels for the central 200 bp. `train.mat` contains 2,200,000 regions, `valid.mat` 4,000 regions, and `test.mat` the chr8/chr9 held-out regions; both forward and reverse-complement examples are represented. The gkm-SVM binary caches are generated from the MAT files by the experiment scripts and do not need separate download.

## 2. hg19 / GRCh37 reference genome

The experiments use UCSC hg19 (GRCh37, February 2009) chromosome FASTAs to extract sequence context and verify reference alleles. The 25 chromosome files (`chr1`–`chr22`, `chrX`, `chrY`, and `chrM`) occupy about 3.16 GB uncompressed. UCSC's compressed files total about 1.0 GB. The current scripts use these root-level paths: `reference_genome/chr1.fa` … `chr22.fa`, `reference_genome/chrX.fa`, and `reference_genome/chrY.fa`.

- Official UCSC directory: <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/>

```bash
mkdir -p reference_genome
for chrom in {1..22} X Y M; do
  curl -fL --retry 3 \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr${chrom}.fa.gz" \
    -o "/tmp/chr${chrom}.fa.gz"
  gzip -dc "/tmp/chr${chrom}.fa.gz" > "reference_genome/chr${chrom}.fa"
done
```

## 3. DeepSEA Supplementary Tables 1–7

These seven files are the paper's supplementary tables and total about 123 MB. Tables 3, 5, and 6 are directly used in current reproductions: Table 3 for DGF/sequence examples, Tables 5–6 for GRASP eQTL and GWAS variants.

- Paper: <https://www.nature.com/articles/nmeth.3547>
- Direct Springer Nature files:

```bash
mkdir -p newdata
base='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects'
for file in \
  41592_2015_BFnmeth3547_MOESM645_ESM.xlsx \
  41592_2015_BFnmeth3547_MOESM646_ESM.xlsx \
  41592_2015_BFnmeth3547_MOESM647_ESM.csv \
  41592_2015_BFnmeth3547_MOESM648_ESM.xlsx \
  41592_2015_BFnmeth3547_MOESM649_ESM.csv \
  41592_2015_BFnmeth3547_MOESM650_ESM.csv \
  41592_2015_BFnmeth3547_MOESM651_ESM.xlsx
do
  curl -fL --retry 3 "$base/$file" -o "newdata/$file"
done
```

## 4. Pretrained model weights

Normal inference downloads and verifies the upstream weights into `~/.cache/deepsea-fuxian3/`; no `models/` directory is needed. The sequence-prediction model is used for chromatin probabilities. The variant-effects model is a separate released checkpoint. Both direct Zenodo endpoints were checked and returned HTTP 200 on 2026-10-07.

- Upstream record: <https://zenodo.org/records/1466993>
- Sequence-prediction model: <https://zenodo.org/records/1466993/files/deepsea_predict.pth?download=1> (MD5 `89e640bf6bdbe1ff165f484d9796efc7`)
- Variant-effects model: <https://zenodo.org/records/1466993/files/deepsea_variant_effects.pth?download=1> (MD5 `35956ab9c28960b5a3693f470fe980c1`)

Optional manual download into the same cache used by `model_assets.py`:

```bash
cache="$HOME/.cache/deepsea-fuxian3"
mkdir -p "$cache"
curl -fL --retry 3 'https://zenodo.org/records/1466993/files/deepsea_predict.pth?download=1' -o "$cache/deepsea_predict.pth"
curl -fL --retry 3 'https://zenodo.org/records/1466993/files/deepsea_variant_effects.pth?download=1' -o "$cache/deepsea_variant_effects.pth"
md5sum "$cache/deepsea_predict.pth" "$cache/deepsea_variant_effects.pth"
```

The locally retrained epoch-53 checkpoint is separate from the pretrained weights. `model_assets.py` fetches it from the project's Hugging Face repository when needed: <https://huggingface.co/aer0vane/reproduce_deepsea>.

## What remains in this checkout

This guide is stored in the `fuxian3/` project root. The data are not deleted from their source; rerunning a section restores files to the root-level paths expected by the code. The training scripts and code remain under `training/` and `experiment/`.
