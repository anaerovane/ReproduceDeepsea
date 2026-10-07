#!/usr/bin/env python3
"""Retry only incomplete Supplementary Fig. 2 TF features.

Features with fewer than 2,000 clean training examples use the maximum
available balanced clean sample. Test subsets use every forward-test positive
and an equal number of randomly selected negatives. A feature with zero test
positives is recorded as not evaluable (no ROC AUC is defined).
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reproduce_suppfig2_all_tf as core

OUT = core.OUT
RETRY_OUT = OUT / "retry_incomplete"
RETRY_INDEX_DIR = RETRY_OUT / "test_indices"
FULL_RESULTS = OUT / "gkm_auc_results.csv"
WORKERS = 4
KERNEL_THREADS = 2


def incomplete_features() -> tuple[int, ...]:
    """Derive retry targets from the recorded full-run statuses."""
    if not FULL_RESULTS.is_file():
        raise FileNotFoundError(f"Full-run result table not found: {FULL_RESULTS}")
    with FULL_RESULTS.open(newline="") as f:
        rows = csv.DictReader(f)
        return tuple(sorted(int(row["feature"]) for row in rows if row.get("status") != "ok"))


def _clean_sample(feature: int, label: int, rng: np.random.Generator) -> list[int]:
    pool = np.flatnonzero(core._TY[:, feature] == label)
    perm = rng.permutation(pool)
    kept: list[int] = []
    cursor = 0
    while len(kept) < core.MAX_TRAIN_PER_CLASS and cursor < len(perm):
        batch = perm[cursor:cursor + 512]
        cursor += len(batch)
        rows = core._TX[batch]
        valid = np.all(rows.sum(axis=2) == 1, axis=1)
        for idx, ok in zip(batch, valid):
            if ok:
                kept.append(int(idx))
                if len(kept) == core.MAX_TRAIN_PER_CLASS:
                    break
    return kept


def _run_one(feature: int, predictor: str) -> dict:
    core.init_worker()

    core.KERNEL_THREADS = KERNEL_THREADS
    seed = core.SEED + feature
    rng = np.random.default_rng(seed)
    train_pos = _clean_sample(feature, 1, rng)
    train_neg = _clean_sample(feature, 0, rng)
    n_train = min(len(train_pos), len(train_neg))
    if n_train == 0:
        return {"feature": feature, "predictor": predictor, "status": "untrainable_no_balanced_clean_samples",
                "n_clean_pos": len(train_pos), "n_clean_neg": len(train_neg)}
    train_pos = rng.choice(train_pos, n_train, replace=False).tolist()
    train_neg = rng.choice(train_neg, n_train, replace=False).tolist()

    y = core._VY[:, feature]
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    if len(positives) == 0:
        return {"feature": feature, "predictor": predictor, "status": "not_evaluable_no_test_positives",
                "n_clean_pos": len(train_pos), "n_clean_neg": len(train_neg),
                "n_test": 0, "n_test_pos": 0}
    n_test_each = min(len(positives), core.MAX_TEST_PER_CLASS)
    test_pos = rng.choice(positives, n_test_each, replace=False)
    test_neg = rng.choice(negatives, n_test_each, replace=False)
    test_idx = np.concatenate([test_pos, test_neg]).astype(np.int32)
    rng.shuffle(test_idx)
    test_y = core._VY[test_idx, feature].astype(np.uint8)

    work = RETRY_OUT / f"feature_{feature}"
    work.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(work / "train_sample_indices.npz",
                        positive=np.asarray(train_pos, dtype=np.int32),
                        negative=np.asarray(train_neg, dtype=np.int32))
    np.save(RETRY_INDEX_DIR / f"feature_{feature}.npy", test_idx)

    seqs: list[str] = []
    for start in range(0, len(train_pos), 128):
        rows = core._TX[train_pos[start:start + 128]]
        seqs.extend("".join(core.BASES[row.argmax(axis=1)]) for row in rows)
    pos_seqs = seqs
    seqs = []
    for start in range(0, len(train_neg), 128):
        rows = core._TX[train_neg[start:start + 128]]
        seqs.extend("".join(core.BASES[row.argmax(axis=1)]) for row in rows)
    neg_seqs = seqs

    records = []
    for start in range(0, len(test_idx), 128):
        idxs = test_idx[start:start + 128]
        rows = core._VX[idxs]
        records.extend((f"test_{int(i)}", "".join(core.BASES[row.argmax(axis=1)]))
                       for i, row in zip(idxs, rows))
    test_fa = work / "test.fa"
    core._write_fasta(test_fa, records)
    aucs = {}
    try:
        for length in (1000, 300):
            model_dir = work / f"gkm{length}"
            sv, alpha = core._train_one(length, pos_seqs, neg_seqs, work, work / "run.log")
            score_file = model_dir / f"test_retry_{len(test_idx)}.scores"
            if not score_file.exists() or score_file.stat().st_size == 0:
                core._run([str(core.BIN_CLASSIFY), "-t", "1", str(test_fa), str(sv), str(alpha), str(score_file)],
                          work / "run.log")
            scores = core._read_scores(score_file)
            expected = {f"test_{int(i)}" for i in test_idx}
            if set(scores) != expected:
                raise RuntimeError(f"score IDs mismatch ({len(scores)}/{len(expected)})")
            vals = np.asarray([scores[f"test_{int(i)}"] for i in test_idx], dtype=np.float64)
            aucs[f"auc_gkm{length}"] = core._auc(test_y, vals)
    finally:
        test_fa.unlink(missing_ok=True)
    return {"feature": feature, "predictor": predictor, "n_train_pos": n_train,
            "n_train_neg": n_train, "n_clean_pos": len(train_pos), "n_clean_neg": len(train_neg),
            "n_test": len(test_idx), "n_test_pos": int(test_y.sum()),
            **aucs, "status": "ok", "error": ""}


def main() -> int:
    features = incomplete_features()
    if not features:
        print("No incomplete features listed in the full-run result table.", flush=True)
        return 0
    RETRY_INDEX_DIR.mkdir(parents=True, exist_ok=True)
    _, all_features = core.prepare_test_cache()
    core.init_worker()
    names = [s.strip() for s in core.PREDICTOR_NAMES.read_text().splitlines() if s.strip()]
    if any(f not in all_features for f in features):
        raise RuntimeError("Retry feature list includes a non-TF feature")
    (RETRY_OUT / "retry_metadata.json").write_text(json.dumps({
        "features": list(features), "workers": WORKERS, "kernel_threads_per_worker": KERNEL_THREADS,
        "training": "maximum available balanced clean examples up to 2000 per class",
        "test": "all available forward-test positives up to 5000 plus equal random negatives",
        "seed": core.SEED, "results_file": "retry_results.csv"
    }, indent=2) + "\n")
    result_file = RETRY_OUT / "retry_results.csv"
    start = time.time()
    with result_file.open("w", newline="") as f:
        fields = ["feature", "predictor", "n_train_pos", "n_train_neg", "n_clean_pos", "n_clean_neg",
                  "n_test", "n_test_pos", "auc_gkm1000", "auc_gkm300", "status", "error"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        f.flush()
        with cf.ProcessPoolExecutor(max_workers=WORKERS) as pool:
            pending = {pool.submit(_run_one, i, names[i]): i for i in features}
            for future in cf.as_completed(pending):
                feature = pending[future]
                try:
                    row = future.result()
                except Exception as e:
                    row = {"feature": feature, "predictor": names[feature], "status": "error", "error": repr(e)}
                writer.writerow(row)
                f.flush()
                print(f"Retry {len(pending)-len([p for p in pending if not p.done()])}/{len(features)} "
                      f"feature={feature} status={row.get('status')} n_train={row.get('n_train_pos')} "
                      f"n_test_pos={row.get('n_test_pos')} auc300={row.get('auc_gkm300')} "
                      f"elapsed_h={(time.time()-start)/3600:.2f}", flush=True)
    core.update_test_index_bundle(all_features)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
