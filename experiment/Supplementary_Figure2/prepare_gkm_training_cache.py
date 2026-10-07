"""Materialize the large byte caches used by the gkm-SVM and Lua workflows.

The upstream bundle supplies compact MAT data. This script converts train.mat
and valid.mat into row-major uint8 binaries expected by the existing runners.
It is only needed after the binary caches have been removed or on a fresh setup.
"""
from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("DEEPSEA_FUXIAN_ROOT", str(HERE.parents[1]))).resolve()
DATA = ROOT / "deepsea_train"
CHUNK_SAMPLES = 64


def materialize_train() -> None:
    source = DATA / "train.mat"
    data_path = DATA / "train_data.bin"
    labels_path = DATA / "train_labels.bin"
    data_tmp = data_path.with_suffix(".bin.partial")
    labels_tmp = labels_path.with_suffix(".bin.partial")

    with h5py.File(source, "r") as mat:
        x = mat["trainxdata"]
        y = mat["traindata"]
        if x.shape[:2] != (1000, 4) or y.shape[0] != 919 or x.shape[2] != y.shape[1]:
            raise ValueError(f"Unexpected train.mat shapes: x={x.shape}, y={y.shape}")
        n = x.shape[2]
        with data_tmp.open("wb") as xf, labels_tmp.open("wb") as yf:
            for start in range(0, n, CHUNK_SAMPLES):
                end = min(start + CHUNK_SAMPLES, n)
                x_chunk = np.asarray(x[:, :, start:end], dtype=np.uint8).transpose(2, 0, 1)
                y_chunk = np.asarray(y[:, start:end], dtype=np.uint8).T
                x_chunk.tofile(xf)
                y_chunk.tofile(yf)
                print(f"train rows {end:,}/{n:,}", flush=True)
    os.replace(data_tmp, data_path)
    os.replace(labels_tmp, labels_path)


def materialize_valid() -> None:
    source = DATA / "valid.mat"
    mat = loadmat(source, variable_names=["validxdata", "validdata"])
    x = np.asarray(mat["validxdata"], dtype=np.uint8)
    y = np.asarray(mat["validdata"], dtype=np.uint8)
    if x.ndim != 3 or x.shape[1:] != (4, 1000) or y.shape != (x.shape[0], 919):
        raise ValueError(f"Unexpected valid.mat shapes: x={x.shape}, y={y.shape}")
    np.ascontiguousarray(x.transpose(0, 2, 1)).tofile(DATA / "valid_data.bin")
    np.ascontiguousarray(y).tofile(DATA / "valid_labels.bin")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    if not (DATA / "train.mat").is_file() or not (DATA / "valid.mat").is_file():
        raise FileNotFoundError(
            f"Expected train.mat and valid.mat in {DATA}; see {DATA / 'DOWNLOAD.md'}"
        )
    materialize_train()
    materialize_valid()
    print("Created train/valid binary caches under deepsea_train/", flush=True)


if __name__ == "__main__":
    main()
