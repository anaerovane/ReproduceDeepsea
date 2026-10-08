"""Recompute Supplementary Figure 7 feature ablations on available cohorts.

Creates separate pretrained and ours figures for GRASP eQTL and GWAS only.
Models are spatial 10-fold OOF XGBoost gblinear classifiers trained on:
  * chromatin effects + four conservation annotations;
  * chromatin effects only;
  * four conservation annotations only.

All metrics are recalculated from fold predictions. Per-fold NPZ files allow
the run to resume. No HGMD or ClinVar proxy is used.
"""
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from huggingface_hub import hf_hub_download

HERE = Path(__file__).resolve().parent
SUPPLEMENTARY_FIGURE6 = HERE.parent / "Supplementary_Figure6"
FIG3 = HERE.parent / "Figure3" / "outputs"
OUT = HERE / "results"
FIGURE_DIR = HERE
GROUPS = ["Random negative SNP", "31kbp negative SNP", "6.3kbp negative SNP",
          "710bp negative SNP", "360bp negative SNP"]
CADD = ["priPhCons", "priPhyloP", "GerpN", "GerpS"]
N_CHROMATIN = 919
MODEL_VERSIONS = ("pretrained", "ours")
HF_REPO_ID = "aer0vane/reproduce_deepsea"


def resolve_artifact(local_path, hub_path):
    """Use a local artifact when present; otherwise fetch the published copy."""
    local_path = Path(local_path)
    if local_path.is_file():
        return local_path
    return Path(hf_hub_download(repo_id=HF_REPO_ID, filename=hub_path,
                               repo_type="model"))


def cohort_path(task):
    name = f"{task}_unsupervised_pretrained.tsv.gz"
    return resolve_artifact(SUPPLEMENTARY_FIGURE6 / name, f"experiment/Supplementary_Figure6/outputs/{name}")


def effect_path(model_name):
    name = f"{model_name}_effects.npy"
    return resolve_artifact(FIG3 / name, f"experiment/Figure3/outputs/{name}")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_cohort(task):


    frame = pd.read_csv(cohort_path(task), sep="\t")
    keep = ["chr", "pos", "ref", "alt", "label", "positive", "fold",
            "distance_bp", "variant_id", *CADD]
    frame = frame[keep].copy()
    for col in CADD:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    n_before = len(frame)
    frame = frame.dropna(subset=CADD).reset_index(drop=True)
    if frame.duplicated(["chr", "pos", "label"]).any():
        raise RuntimeError(f"{task}: duplicate position/label records in source cohort")
    print(f"{task}: common complete cases {len(frame):,}/{n_before:,}; "
          f"positive={int(frame.positive.sum()):,}; "
          f"random-negative={int(frame.label.eq(GROUPS[0]).sum()):,}", flush=True)
    return frame


def feature_matrix(frame, model_name, kind):
    ids = frame.variant_id.to_numpy(dtype=np.int64)
    if kind == "conservation":
        return frame[CADD].to_numpy(dtype=np.float32)
    effects = np.load(effect_path(model_name), mmap_mode="r")
    if effects.ndim != 2 or effects.shape[1] != 2 * N_CHROMATIN:
        raise RuntimeError(f"Unexpected {model_name} effect shape: {effects.shape}")
    raw = np.asarray(effects[ids], dtype=np.float32)
    chromatin = np.abs(raw[:, :N_CHROMATIN]) * np.abs(raw[:, N_CHROMATIN:])
    if kind == "chromatin":
        return chromatin
    if kind == "combined":
        return np.concatenate([chromatin, frame[CADD].to_numpy(dtype=np.float32)], axis=1)
    raise ValueError(kind)


def fit_predict_fold(x, frame, fold, task, model_tag, model_version, out):
    cache_tag = model_tag if model_version == "shared" else f"{model_version}_{model_tag}"
    cache = out / f"{task}_{cache_tag}_fold{fold:02d}.npz"
    test_mask = frame.fold.to_numpy() == fold
    test_ix = np.flatnonzero(test_mask)
    if cache.exists():
        z = np.load(cache)
        if np.array_equal(z["test_indices"], test_ix) and len(z["probabilities"]) == len(test_ix):
            print(f"{task}/{cache_tag}: restored fold {fold}", flush=True)
            return test_ix, z["probabilities"]
        raise RuntimeError(f"Stale fold cache: {cache}")

    train_base = frame.positive.to_numpy() | frame.label.eq(GROUPS[0]).to_numpy()
    train_mask = (frame.fold.to_numpy() != fold) & train_base
    train_ix = np.flatnonzero(train_mask)
    y = frame.positive.to_numpy(dtype=np.int8)[train_ix]
    weights = np.where(y == 1, len(y) / (2 * max((y == 1).sum(), 1)),
                       len(y) / (2 * max((y == 0).sum(), 1))).astype(np.float32)
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x[train_ix]).astype(np.float32, copy=False)
    x_test = scaler.transform(x[test_ix]).astype(np.float32, copy=False)
    penalty_scale = float(weights.sum())

    class Progress(xgb.callback.TrainingCallback):
        def after_iteration(self, model, epoch, evals_log):
            if (epoch + 1) % 20 == 0 or epoch + 1 == 100:
                print(f"{task}/{cache_tag}/fold{fold}: {epoch + 1}/100", flush=True)
            return False

    classifier = xgb.XGBClassifier(
        booster="gblinear", updater="coord_descent", feature_selector="cyclic",
        device="cuda" if torch.cuda.is_available() else "cpu",
        reg_alpha=0.0, reg_lambda=10.0 / penalty_scale,
        n_estimators=100, learning_rate=0.1, callbacks=[Progress()],
        n_jobs=8, random_state=42, eval_metric="logloss",
    )
    classifier.fit(x_train, y, sample_weight=weights)
    probabilities = classifier.predict_proba(x_test)[:, 1].astype(np.float32)
    tmp = cache.with_suffix(".npz.tmp")
    with open(tmp, "wb") as f:
        np.savez_compressed(f, test_indices=test_ix, probabilities=probabilities,
                            train_n=np.array([len(train_ix)]), test_n=np.array([len(test_ix)]))
    os.replace(tmp, cache)
    print(f"{task}/{cache_tag}: completed fold {fold}; "
          f"train={len(train_ix):,}, test={len(test_ix):,}", flush=True)
    del classifier, x_train, x_test
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return test_ix, probabilities


def evaluate(frame, x, task, model_tag, out, model_version="shared"):
    pred = np.full(len(frame), np.nan, dtype=np.float32)
    for fold in range(10):
        test_ix, p = fit_predict_fold(x, frame, fold, task, model_tag,
                                      model_version, out)
        pred[test_ix] = p
    if not np.isfinite(pred).all():
        missing = int((~np.isfinite(pred)).sum())
        raise RuntimeError(f"{task}/{model_version}/{model_tag}: {missing} rows lack OOF predictions")
    rows = []
    y = frame.positive.to_numpy(dtype=np.int8)
    for group in GROUPS:
        ix = frame.positive.to_numpy() | frame.label.eq(group).to_numpy()
        neg = frame.label.eq(group).to_numpy()
        rows.append({
            "task": task, "feature_set": model_tag, "negative_group": group,
            "mean_distance_bp": float(frame.loc[neg, "distance_bp"].mean()),
            "n_positive": int((ix & (y == 1)).sum()),
            "n_negative": int(neg.sum()),
            "auc": float(roc_auc_score(y[ix], pred[ix])),
        })
    pred_frame = frame[["chr", "pos", "ref", "alt", "label", "positive",
                        "fold", "distance_bp", "variant_id"]].copy()
    pred_frame["oof_probability"] = pred
    oof_tag = model_tag if model_version == "shared" else f"{model_version}_{model_tag}"
    pred_frame.to_csv(out / f"{task}_{oof_tag}_oof.tsv.gz", sep="\t", index=False)
    return rows


def make_plots(metrics, out):
    task_titles = {"eqtl": "GRASP eQTL (noncoding)",
                   "gwas": "GWAS Catalog (noncoding)"}
    styles = {
        "combined": ("#13a9d2", "-", "DeepSEA"),
        "chromatin": ("#13a9d2", "--", "DeepSEA (Predicted chromatin effect only)"),
        "conservation": ("#13a9d2", ":", "DeepSEA (Evolutionary information only)"),
    }
    for model_name in MODEL_VERSIONS:
        fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.65),
                                 gridspec_kw={"wspace": 0.28})
        display_name = "Pretrained" if model_name == "pretrained" else "Ours"
        for ax, task in zip(axes, ("eqtl", "gwas")):
            task_rows = metrics[(metrics.task == task) &
                                (metrics.model_version == model_name)]
            for feature_set in ("combined", "chromatin", "conservation"):
                subset = task_rows[task_rows.feature_set == feature_set].set_index("negative_group")
                subset = subset.reindex(GROUPS)
                color, linestyle, label = styles[feature_set]
                ax.plot(np.arange(len(GROUPS)), subset.auc.to_numpy(dtype=float),
                        linestyle=linestyle, color=color, lw=1.5, label=label)
            base = task_rows[task_rows.feature_set == "combined"].set_index("negative_group").reindex(GROUPS)
            labels = ["All"] + [f"{x:,.0f}" for x in base.mean_distance_bp.to_numpy()[1:]]
            n_pos = int(base.n_positive.iloc[0])
            ax.set_title(f"{task_titles[task]}\n($n$ = {n_pos:,})", fontsize=9, pad=4)
            ax.set_xticks(np.arange(len(GROUPS)), labels)
            ax.set_xlabel("Negative SNP group (mean distance, bp)", fontsize=8, labelpad=3)
            ax.set_ylabel("AUC", fontsize=8)
            ax.set_ylim(0.5, 0.76)
            ax.tick_params(axis="both", labelsize=7.5, length=2, pad=2)
            ax.spines[["top", "right"]].set_visible(False)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, fontsize=7.3,
                   loc="center left", bbox_to_anchor=(0.79, 0.52))
        fig.suptitle(display_name, fontsize=10, fontweight="semibold", y=0.99)
        fig.subplots_adjust(left=0.075, right=0.77, bottom=0.22, top=0.80)
        for ext in ("png", "pdf"):
            fig.savefig(FIGURE_DIR / f"suppfig7_{model_name}.{ext}", dpi=300,
                        bbox_inches="tight")
        plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(8)
    torch.backends.cudnn.benchmark = True
    manifest_path = OUT / "run_manifest.json"
    input_files = [cohort_path(task) for task in ("eqtl", "gwas")]
    input_files += [effect_path(name) for name in MODEL_VERSIONS]
    manifest = {
        "description": "Supplementary Figure 7 feature ablation; GRASP and GWAS only",
        "feature_sets": {
            "combined": "919 DeepSEA effect magnitudes + CADD GRCh37-v1.4 inclAnno priPhCons, priPhyloP, GerpN, GerpS",
            "chromatin": "919 DeepSEA effect magnitudes",
            "conservation": "four CADD GRCh37-v1.4 inclAnno conservation annotations",
        },
        "model": "XGBoost gblinear; weighted binary logistic; 10 supplied spatial folds; per-fold StandardScaler; n_estimators=100; reg_lambda=10",
        "negative_groups": GROUPS,
        "missing_data": "retain common rows complete for all four CADD fields; all random-negative rows retained",
        "excluded": ["HGMD (source cohort unavailable)", "ClinVar proxy"],
        "inputs_sha256": {p.name: sha256(p) for p in input_files},
        "code_sha256": sha256(Path(__file__)),
        "cache_layout": "v2; DeepSEA model version included in chromatin/combined fold cache keys",
        "device": "cuda" if torch.cuda.is_available() else "cpu",
    }
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        stable_keys = ["description", "feature_sets", "model", "negative_groups",
                       "missing_data", "excluded", "inputs_sha256", "device"]
        if any(previous.get(key) != manifest.get(key) for key in stable_keys):
            raise RuntimeError("Existing output cache has different inputs/config; use a new folder")
        manifest["started_at"] = previous.get("started_at", time.strftime("%Y-%m-%dT%H:%M:%S%z"))
        if previous.get("code_sha256") != manifest["code_sha256"]:
            manifest["previous_code_sha256"] = previous.get("code_sha256")
            manifest["cache_correction"] = "Pretrained and ours caches are now separated; prior unlabeled DeepSEA fold files are treated as pretrained only."
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        manifest = previous
    else:
        manifest["started_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    all_rows = []


    conservation_by_task = {}
    for task in ("eqtl", "gwas"):
        frame = read_cohort(task)
        x = feature_matrix(frame, "pretrained", "conservation")
        conservation_by_task[task] = evaluate(frame, x, task, "conservation", OUT)
        del x

    for model_name in MODEL_VERSIONS:
        for task in ("eqtl", "gwas"):
            frame = read_cohort(task)
            for feature_set in ("chromatin", "combined"):
                print(f"Starting {model_name}/{task}/{feature_set}", flush=True)
                x = feature_matrix(frame, model_name, feature_set)
                rows = evaluate(frame, x, task, feature_set, OUT,
                                model_version=model_name)
                for row in rows:
                    row["model_version"] = model_name
                all_rows.extend(rows)
                del x
            for row in conservation_by_task[task]:
                row = dict(row)
                row["model_version"] = model_name
                all_rows.append(row)
    metrics = pd.DataFrame(all_rows)
    metrics.to_csv(OUT / "auc_by_group.csv", index=False)
    make_plots(metrics, OUT)
    manifest["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    manifest["auc_sha256"] = sha256(OUT / "auc_by_group.csv")
    manifest["figures_sha256"] = {
        p.name: sha256(p) for p in sorted(FIGURE_DIR.glob("suppfig7_*"))
        if p.suffix in {".png", ".pdf"}
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(metrics.to_string(index=False), flush=True)
    print(f"Wrote figures and metrics under {OUT}", flush=True)


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        metrics = pd.read_csv(OUT / "auc_by_group.csv")
        make_plots(metrics, OUT)
        for model_name in MODEL_VERSIONS:
            for ext in ("png", "pdf"):
                name = f"suppfig7_{model_name}.{ext}"
                shutil.copy2(HERE / name, HERE.parent / name)
        manifest_path = OUT / "run_manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            manifest["figure_model_labels"] = {"pretrained": "Pretrained", "ours": "Ours"}
            manifest["figures_sha256"] = {
                p.name: sha256(p) for p in sorted(HERE.glob("suppfig7_*"))
                if p.suffix in {".png", ".pdf"}
            }
            manifest["figures_location"] = "experiment/Supplementary_Figure7 root"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print("Redrew labeled figures from existing AUC table; no model evaluation was run.")
    else:
        main()
