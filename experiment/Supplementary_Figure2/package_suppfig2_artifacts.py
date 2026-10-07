#!/usr/bin/env python3
"""Merge per-feature outputs, publish the archive to HF, and leave no local copy."""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from huggingface_hub import HfApi, hf_hub_download

HERE = Path(__file__).resolve().parent
OUT = HERE / "suppfig2_tf_fulltrain"
ARCHIVE_NAME = "feature_models_and_indices.tar.zst"
ARCHIVE = OUT / ARCHIVE_NAME
HF_REPO_ID = os.environ.get("DEEPSEA_HF_REPO", "aer0vane/reproduce_deepsea")
HF_ARCHIVE_PATH = f"experiment/Supplementary_Figure2/{ARCHIVE_NAME}"
PACKED_DIRS = ("retry_incomplete", "test_indices")


def run(args: list[str]) -> str:
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch_previous_archive(destination: Path) -> Path:
    if ARCHIVE.is_file():
        return ARCHIVE
    path = hf_hub_download(
        repo_id=HF_REPO_ID,
        filename=HF_ARCHIVE_PATH,
        repo_type="model",
        local_dir=destination / "download",
    )
    return Path(path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pending = [p for p in OUT.glob("feature_*") if p.is_dir()]
    pending += [OUT / name for name in PACKED_DIRS if (OUT / name).is_dir()]
    if not pending:
        print("No expanded run folders to package; no large local archive created.")
        return

    api = HfApi()
    with tempfile.TemporaryDirectory(prefix="suppfig2-package-") as temp:
        temp_root = Path(temp)
        stage = temp_root / "stage"
        stage.mkdir()

        previous = fetch_previous_archive(temp_root)
        listing = run(["tar", "-I", "zstd", "-tf", str(previous)]).splitlines()
        strip = 1 if listing and all(x.startswith("suppfig2_tf_fulltrain/") for x in listing) else 0
        extract = ["tar", "-I", "zstd", "-xf", str(previous), "-C", str(stage)]
        if strip:
            extract += ["--strip-components=1"]
        subprocess.run(extract, check=True)

        # Merge new run outputs over the archived files while keeping previous TFs.
        for source in pending:
            shutil.copytree(source, stage / source.name, dirs_exist_ok=True)
        names = sorted(p.name for p in stage.glob("feature_*") if p.is_dir())
        names += [name for name in PACKED_DIRS if (stage / name).is_dir()]

        # Retry-specific test indices take precedence in the directly usable bundle.
        index_dir = stage / "test_indices"
        retry_dir = stage / "retry_incomplete" / "test_indices"
        arrays = {}
        for path in sorted(index_dir.glob("feature_*.npy")) if index_dir.exists() else []:
            retry_path = retry_dir / path.name
            arrays[path.stem] = np.load(
                retry_path if retry_path.exists() else path, allow_pickle=False
            )
        if arrays:
            np.savez_compressed(OUT / "test_indices.npz", **arrays)

        if not names:
            print("No model/index artifacts found; expanded folders retained.")
            return

        packed = temp_root / ARCHIVE_NAME
        subprocess.run([
            "tar", "--use-compress-program=zstd -T0 -3", "-cf", str(packed),
            "-C", str(stage), "--", *names,
        ], check=True)
        subprocess.run(["zstd", "-t", str(packed)], check=True, stdout=subprocess.DEVNULL)
        run(["tar", "-I", "zstd", "-tf", str(packed)])

        sha256 = digest(packed)
        api.upload_file(
            path_or_fileobj=str(packed),
            path_in_repo=HF_ARCHIVE_PATH,
            repo_id=HF_REPO_ID,
            repo_type="model",
            commit_message="Update Supplementary Figure 2 gkm-SVM model archive",
        )
        info = api.model_info(HF_REPO_ID, files_metadata=True)
        remote = next((x for x in info.siblings if x.rfilename == HF_ARCHIVE_PATH), None)
        remote_hash = getattr(getattr(remote, "lfs", None), "sha256", None)
        if remote is None or remote.size != packed.stat().st_size or remote_hash != sha256:
            raise RuntimeError("Hugging Face archive verification failed; expanded folders retained")

        # Remove extracted outputs only after the matching remote file is confirmed.
        for path in pending:
            if path.is_dir():
                shutil.rmtree(path)
        ARCHIVE.unlink(missing_ok=True)
        print(f"Uploaded and verified {ARCHIVE_NAME} on {HF_REPO_ID}; local archive removed.")


if __name__ == "__main__":
    main()
