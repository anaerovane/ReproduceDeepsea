#!/usr/bin/env python3
"""Package per-feature run folders into a flat, lossless zstd tar archive."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "suppfig2_tf_fulltrain"
ARCHIVE = OUT / "feature_models_and_indices.tar.zst"
PACKED_DIRS = ("retry_incomplete", "test_indices")


def run(args: list[str]) -> str:
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="suppfig2-package-") as temp:
        stage = Path(temp)

        # Merge the previous archive so packaging a retry never drops earlier models.
        if ARCHIVE.exists():
            listing = run(["tar", "-I", "zstd", "-tf", str(ARCHIVE)]).splitlines()
            strip = 1 if listing and all(x.startswith("suppfig2_tf_fulltrain/") for x in listing) else 0
            command = ["tar", "-I", "zstd", "-xf", str(ARCHIVE), "-C", str(stage)]
            if strip:
                command += ["--strip-components=1"]
            subprocess.run(command, check=True)

        names = sorted(p.name for p in stage.glob("feature_*") if p.is_dir())
        names += [n for n in PACKED_DIRS if (stage / n).is_dir()]
        for source in sorted(OUT.glob("feature_*")):
            if source.is_dir():
                names.append(source.name)
                shutil.copytree(source, stage / source.name, dirs_exist_ok=True)
        for name in PACKED_DIRS:
            source = OUT / name
            if source.is_dir():
                if name not in names:
                    names.append(name)
                shutil.copytree(source, stage / name, dirs_exist_ok=True)
        names = sorted(set(names))

        # Preserve retry indices in the flat bundle, overriding original indices.
        index_dir = stage / "test_indices"
        retry_dir = stage / "retry_incomplete" / "test_indices"
        arrays = {}
        for p in sorted(index_dir.glob("feature_*.npy")) if index_dir.exists() else []:
            feature = p.stem
            retry = retry_dir / f"{feature}.npy"
            arrays[feature] = np.load(retry if retry.exists() else p, allow_pickle=False)
        if arrays:
            np.savez_compressed(OUT / "test_indices.npz", **arrays)

        if not names:
            print("No per-feature folders found; existing archive left unchanged.")
            return
        temp_archive = OUT / ".feature_models_and_indices.tar.zst.tmp"
        temp_archive.unlink(missing_ok=True)
        subprocess.run([
            "tar", "--use-compress-program=zstd -T0 -3", "-cf", str(temp_archive),
            "-C", str(stage), "--", *names,
        ], check=True)
        subprocess.run(["zstd", "-t", str(temp_archive)], check=True)
        run(["tar", "-I", "zstd", "-tf", str(temp_archive)])
        temp_archive.replace(ARCHIVE)

        for p in list(OUT.glob("feature_*")):
            if p.is_dir():
                shutil.rmtree(p)
        for name in PACKED_DIRS:
            p = OUT / name
            if p.is_dir():
                shutil.rmtree(p)
        print(f"Packaged {len(names)} top-level entries into {ARCHIVE.name}")


if __name__ == "__main__":
    main()
