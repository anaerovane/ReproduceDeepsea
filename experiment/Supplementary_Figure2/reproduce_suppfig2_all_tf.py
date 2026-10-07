#!/usr/bin/env python3
"""Reproduce the gkm-SVM TF comparison for DeepSEA Supplementary Fig. 2.

Uses the true predictor order, full train split for sampling, both 1000bp and
center-300bp gkm-SVM models, and correctly serialized forward test.mat rows.
Per-feature balanced test subsets (up to 5,000 positives plus the same number
of negatives) keep the run tractable; the subset indices are saved verbatim.
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("DEEPSEA_FUXIAN_ROOT", str(HERE.parents[1]))).resolve()
DATA = ROOT / "deepsea_train"
OUT = HERE / "results"
BIN_KERNEL = HERE / "gkmsvm_kernel"
BIN_TRAIN = HERE / "gkmsvm_train"
BIN_CLASSIFY = HERE / "gkmsvm_classify"
N_FEATURES = 919
N_TRAIN = 4_400_000
WIDTH = 1000
MAX_TRAIN_PER_CLASS = int(os.environ.get("GKM_TRAIN_PER_CLASS", "2000"))
MAX_TEST_PER_CLASS = int(os.environ.get("GKM_TEST_PER_CLASS", "5000"))
KERNEL_THREADS = int(os.environ.get("GKM_KERNEL_THREADS", "4"))
MIN_TEST_POS = 10
SEED = 20261006
N_WORKERS = int(os.environ.get("GKM_WORKERS", "8"))
BASES = np.array(list("AGCT"))

TRAIN_X = DATA / "train_data.bin"
TRAIN_Y = DATA / "train_labels.bin"
TEST_X = OUT / "test_forward_agct_1000x4.bin"
TEST_Y = OUT / "test_forward_labels_rowmajor.bin"
TEST_INDEX_DIR = OUT / "test_indices"
TEST_INDEX_BUNDLE = OUT / "test_indices.npz"
PREDICTOR_NAMES = ROOT / "predictor_names.txt"

_TX = None
_TY = None
_VX = None
_VY = None


def predictor_tf_indices() -> tuple[list[str], list[int]]:
    names = [s.strip() for s in PREDICTOR_NAMES.read_text().splitlines() if s.strip()]
    if len(names) != N_FEATURES:
        raise RuntimeError(f"Expected {N_FEATURES} predictor names, got {len(names)}")
    histone = re.compile(r"_(?:H2[A-Za-z0-9]+|H3[A-Za-z0-9]+|H4[A-Za-z0-9]+)_")
    indices = [i for i, name in enumerate(names)
               if "_DNase_" not in name and not histone.search("_" + name + "_")]
    if len(indices) != 690:
        raise RuntimeError(f"Expected 690 TF predictors, found {len(indices)}")
    return names, indices


def prepare_test_cache() -> tuple[list[str], list[int]]:
    names, tf_indices = predictor_tf_indices()
    TEST_INDEX_DIR.mkdir(parents=True, exist_ok=True)
    if not TEST_X.exists() or not TEST_Y.exists():
        from scipy.io import loadmat
        mat = loadmat(DATA / "test.mat", variable_names=["testxdata", "testdata"])
        x = mat["testxdata"]
        y = mat["testdata"]
        if x.shape != (455024, 4, WIDTH) or y.shape != (455024, N_FEATURES):
            raise RuntimeError(f"Unexpected test.mat dimensions: {x.shape}, {y.shape}")
        n_forward = x.shape[0] // 2
        # DeepSEA test.mat is AGCT, channels first. Store row-major (sample, position, channel).
        x[:n_forward].transpose(0, 2, 1).tofile(TEST_X)
        y[:n_forward].tofile(TEST_Y)
        del mat, x, y

    n_test = 227_512
    vx = np.memmap(TEST_X, dtype=np.uint8, mode="r", shape=(n_test, WIDTH, 4))
    vy = np.memmap(TEST_Y, dtype=np.uint8, mode="r", shape=(n_test, N_FEATURES))
    # Strong serialization guard: source row 0 reconstructed from the raw cache.
    from scipy.io import loadmat
    m = loadmat(DATA / "test.mat", variable_names=["testxdata", "testdata"])
    source_seq = "".join(BASES[np.argmax(m["testxdata"][0], axis=0)])
    cache_seq = "".join(BASES[np.argmax(vx[0], axis=1)])
    if source_seq != cache_seq or not np.array_equal(vy[0], m["testdata"][0]):
        raise RuntimeError("Correctly serialized test cache failed row-0 round-trip")
    del m

    # Persist deterministic, feature-stratified held-out indices from the forward half.
    rng = np.random.default_rng(SEED)
    for feature in tf_indices:
        path = TEST_INDEX_DIR / f"feature_{feature}.npy"
        if path.exists():
            continue
        positives = np.flatnonzero(vy[:, feature] == 1)
        negatives = np.flatnonzero(vy[:, feature] == 0)
        if len(positives) < MIN_TEST_POS:
            np.save(path, np.empty(0, dtype=np.int32))
            continue
        n = min(len(positives), MAX_TEST_PER_CLASS)
        p = rng.choice(positives, n, replace=False)
        q = rng.choice(negatives, n, replace=False)
        idx = np.concatenate((p, q)).astype(np.int32)
        rng.shuffle(idx)
        np.save(path, idx)
    update_test_index_bundle(tf_indices)
    return names, tf_indices


def update_test_index_bundle(features: list[int]) -> None:
    """Write a flat index bundle, preferring retry-specific held-out rows."""
    retry_dir = OUT / "retry_incomplete" / "test_indices"
    arrays = {}
    for feature in features:
        retry_path = retry_dir / f"feature_{feature}.npy"
        path = retry_path if retry_path.exists() else TEST_INDEX_DIR / f"feature_{feature}.npy"
        if path.exists():
            arrays[f"feature_{feature}"] = np.load(path, allow_pickle=False).astype(np.int32)
    np.savez_compressed(TEST_INDEX_BUNDLE, **arrays)


def init_worker() -> None:
    global _TX, _TY, _VX, _VY
    _TX = np.memmap(TRAIN_X, dtype=np.uint8, mode="r", shape=(N_TRAIN, WIDTH, 4))
    _TY = np.memmap(TRAIN_Y, dtype=np.uint8, mode="r", shape=(N_TRAIN, N_FEATURES))
    _VX = np.memmap(TEST_X, dtype=np.uint8, mode="r", shape=(227_512, WIDTH, 4))
    _VY = np.memmap(TEST_Y, dtype=np.uint8, mode="r", shape=(227_512, N_FEATURES))


def _write_fasta(path: Path, records) -> None:
    with path.open("w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def _read_scores(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    with path.open() as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) == 2:
                result[p[0]] = float(p[1])
    return result


def _auc(y: np.ndarray, s: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    if len(y) != len(s) or len(np.unique(y)) != 2:
        return float("nan")
    return float(roc_auc_score(y, s))


def _run(cmd: list[str], log: Path, timeout: int = 7200) -> None:
    with log.open("a") as f:
        f.write("\n$ " + " ".join(cmd) + "\n")
        f.flush()
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd[0]} (see {log})")


def _sample_training(feature: int, work: Path) -> tuple[list[str], list[str], list[int], list[int]]:
    assert _TX is not None and _TY is not None
    rng = np.random.default_rng(SEED + feature)
    classes: list[list[int]] = []
    for label in (1, 0):
        pool = np.flatnonzero(_TY[:, feature] == label)
        perm = rng.permutation(pool)
        kept: list[int] = []
        cursor = 0
        while len(kept) < MAX_TRAIN_PER_CLASS and cursor < len(perm):
            batch = perm[cursor:cursor + 512]
            cursor += len(batch)
            rows = _TX[batch]
            valid = np.all(rows.sum(axis=2) == 1, axis=1)
            for idx, ok in zip(batch, valid):
                if ok:
                    kept.append(int(idx))
                    if len(kept) == MAX_TRAIN_PER_CLASS:
                        break
        if len(kept) < MAX_TRAIN_PER_CLASS:
            raise RuntimeError(f"feature {feature} label {label}: only {len(kept)} clean examples")
        classes.append(kept)

    pos_idx, neg_idx = classes
    pos_seqs: list[str] = []
    neg_seqs: list[str] = []
    for idxs, dest in ((pos_idx, pos_seqs), (neg_idx, neg_seqs)):
        for start in range(0, len(idxs), 128):
            rows = _TX[idxs[start:start + 128]]
            dest.extend("".join(BASES[row.argmax(axis=1)]) for row in rows)
    np.savez_compressed(work / "train_sample_indices.npz",
                        positive=np.asarray(pos_idx, dtype=np.int32),
                        negative=np.asarray(neg_idx, dtype=np.int32))
    return pos_seqs, neg_seqs, pos_idx, neg_idx


def _train_one(length: int, pos: list[str], neg: list[str], work: Path, log: Path) -> tuple[Path, Path]:
    dest = work / f"gkm{length}"
    dest.mkdir(exist_ok=True)
    sv = dest / "model_svseq.fa"
    alpha = dest / "model_svalpha.out"
    if sv.exists() and alpha.exists() and sv.stat().st_size and alpha.stat().st_size:
        return sv, alpha
    pfa, nfa = work / f"train_pos_{length}.fa", work / f"train_neg_{length}.fa"
    if length == WIDTH:
        pseq, nseq = pos, neg
    else:
        lo = (WIDTH - length) // 2
        hi = lo + length
        pseq = [s[lo:hi] for s in pos]
        nseq = [s[lo:hi] for s in neg]
    _write_fasta(pfa, ((f"p{i}", s) for i, s in enumerate(pseq)))
    _write_fasta(nfa, ((f"n{i}", s) for i, s in enumerate(nseq)))
    kernel = dest / "kernel.txt"
    prefix = dest / "model"
    _run([str(BIN_KERNEL), "-l", "10", "-k", "6", "-d", "3", "-T", str(KERNEL_THREADS), str(pfa), str(nfa), str(kernel)], log)
    _run([str(BIN_TRAIN), str(kernel), str(pfa), str(nfa), str(prefix)], log)
    kernel.unlink(missing_ok=True)
    pfa.unlink(missing_ok=True)
    nfa.unlink(missing_ok=True)
    if not sv.exists() or not alpha.exists():
        raise RuntimeError(f"gkm-SVM did not produce model files for {dest}")
    return sv, alpha


def run_feature(feature: int) -> dict:
    assert _VX is not None and _VY is not None
    work = OUT / f"feature_{feature}"
    work.mkdir(parents=True, exist_ok=True)
    log = work / "run.log"
    sample_idx_file = TEST_INDEX_DIR / f"feature_{feature}.npy"
    test_idx = np.load(sample_idx_file)
    if len(test_idx) == 0:
        return {"feature": feature, "status": "too_few_test_positives"}
    test_y = _VY[test_idx, feature].astype(np.uint8)
    if len(np.unique(test_y)) != 2:
        return {"feature": feature, "status": "single_test_class"}

    if (work / "train_sample_indices.npz").exists():
        saved = np.load(work / "train_sample_indices.npz")
        pos_idx, neg_idx = saved["positive"].tolist(), saved["negative"].tolist()
        pos_seqs, neg_seqs = [], []
        for idxs, dest in ((pos_idx, pos_seqs), (neg_idx, neg_seqs)):
            for start in range(0, len(idxs), 128):
                rows = _TX[idxs[start:start + 128]]
                dest.extend("".join(BASES[row.argmax(axis=1)]) for row in rows)
    else:
        pos_seqs, neg_seqs, pos_idx, neg_idx = _sample_training(feature, work)

    result = {"feature": feature, "n_train_pos": len(pos_idx), "n_train_neg": len(neg_idx),
              "n_test": len(test_idx), "n_test_pos": int(test_y.sum()), "status": "ok"}
    test_records = []
    for start in range(0, len(test_idx), 128):
        idxs = test_idx[start:start + 128]
        rows = _VX[idxs]
        test_records.extend((f"test_{int(i)}", "".join(BASES[row.argmax(axis=1)]))
                            for i, row in zip(idxs, rows))
    full_fa = work / "test_1000.fa"
    short_fa = work / "test_300.fa"
    _write_fasta(full_fa, test_records)
    _write_fasta(short_fa, ((name, seq[(WIDTH-300)//2:(WIDTH+300)//2]) for name, seq in test_records))

    try:
        for length, test_fa in ((1000, full_fa), (300, short_fa)):
            sv, alpha = _train_one(length, pos_seqs, neg_seqs, work, log)
            scores_file = work / f"gkm{length}" / f"test_{len(test_idx)}.scores"
            if not scores_file.exists() or scores_file.stat().st_size == 0:
                _run([str(BIN_CLASSIFY), "-t", "1", str(test_fa), str(sv), str(alpha), str(scores_file)], log)
            scores = _read_scores(scores_file)
            expected = {f"test_{int(i)}" for i in test_idx}
            if set(scores) != expected:
                raise RuntimeError(f"feature {feature}, {length}bp score IDs mismatch: {len(scores)}/{len(expected)}")
            vals = np.asarray([scores[f"test_{int(i)}"] for i in test_idx], dtype=np.float64)
            auc = _auc(test_y, vals)
            result[f"auc_gkm{length}"] = auc
        return result
    finally:
        full_fa.unlink(missing_ok=True)
        short_fa.unlink(missing_ok=True)


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", nargs="*", type=int, help="optional predictor indices")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if not TRAIN_X.is_file() or not TRAIN_Y.is_file():
        print("Expanded training cache is missing; rebuilding it from train.mat (about 20 GB).", flush=True)
        subprocess.run(
            [sys.executable, str(HERE / "prepare_gkm_training_cache.py")],
            cwd=ROOT,
            env={**os.environ, "DEEPSEA_FUXIAN_ROOT": str(ROOT)},
            check=True,
        )
    names, features = prepare_test_cache()
    if args.features is not None:
        bad = sorted(set(args.features) - set(features))
        if bad:
            raise ValueError(f"Requested non-TF feature indices: {bad}")
        features = list(args.features)
    (OUT / "run_metadata.json").write_text(json.dumps({
        "objective": "Supplementary Figure 2 gkm-SVM comparison, TF features",
        "n_tf_features": len(features), "tf_feature_indices": [features[0], features[-1]],
        "first_tf": names[features[0]], "last_tf": names[features[-1]],
        "train_draw_per_class": MAX_TRAIN_PER_CLASS,
        "test_sample": "per TF, all forward-test positives up to 5000 and an equal random negative sample",
        "test_rows": 227512, "random_seed": SEED, "workers": N_WORKERS,
        "alphabet": "AGCT", "models": ["1000bp", "center 300bp"]
    }, indent=2) + "\n")
    print(f"TF features {len(features)}: {features[0]}..{features[-1]} ({names[features[0]]} .. {names[features[-1]]})", flush=True)
    print(f"Balanced held-out test indices prepared; workers={N_WORKERS}", flush=True)
    results_path = OUT / "gkm_auc_results.csv"
    done: set[int] = set()
    if results_path.exists():
        with results_path.open() as f:
            done = {int(row["feature"]) for row in csv.DictReader(f) if row.get("status") == "ok"}
    todo = [i for i in features if i not in done]
    print(f"Already complete: {len(done)}/{len(features)}; remaining {len(todo)}", flush=True)
    if not todo:
        return 0
    exists = results_path.exists()
    with results_path.open("a", newline="") as out_csv:
        fields = ["feature", "predictor", "n_train_pos", "n_train_neg", "n_test", "n_test_pos",
                  "auc_gkm1000", "auc_gkm300", "status", "error"]
        writer = csv.DictWriter(out_csv, fieldnames=fields)
        if not exists:
            writer.writeheader()
            out_csv.flush()
        start = time.time()
        completed = 0
        with cf.ProcessPoolExecutor(max_workers=N_WORKERS, initializer=init_worker) as pool:
            pending = iter(todo)
            futures = {pool.submit(run_feature, i): i for i in
                       (next(pending, None) for _ in range(min(N_WORKERS, len(todo))))}
            futures = {f: i for f, i in futures.items() if i is not None}
            while futures:
                fut = next(cf.as_completed(futures))
                i = futures.pop(fut)
                try:
                    row = fut.result()
                    row["predictor"] = names[i]
                    row["error"] = ""
                except Exception as e:
                    row = {"feature": i, "predictor": names[i], "status": "error", "error": repr(e)}
                writer.writerow(row)
                out_csv.flush()
                completed += 1
                elapsed = time.time() - start
                rate = completed / elapsed if elapsed else 0
                eta = (len(todo) - completed) / rate if rate else float("inf")
                print(f"Progress {len(done)+completed}/{len(features)} | feature {i} {names[i]} | "
                      f"status={row.get('status')} AUC300={row.get('auc_gkm300')} "
                      f"elapsed={elapsed/3600:.2f}h ETA={eta/3600:.1f}h", flush=True)
                nxt = next(pending, None)
                if nxt is not None:
                    futures[pool.submit(run_feature, nxt)] = nxt
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
