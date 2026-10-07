#!/usr/bin/env python3
"""Download the public GEO count tracks and UCSC hg18→hg19 chain on demand."""
from __future__ import annotations

import csv
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
MANIFEST = INPUTS / "GEO_GSE47991_download_manifest.csv"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"


def download(url: str, target: Path, expected_bytes: int | None = None) -> None:
    if target.exists() and (expected_bytes is None or target.stat().st_size == expected_bytes):
        print(f"present: {target.name}")
        return
    temp = target.with_suffix(target.suffix + ".partial")
    total = 0
    with urlopen(url, timeout=90) as response, temp.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
            total += len(chunk)
    if expected_bytes is not None and total != expected_bytes:
        temp.unlink(missing_ok=True)
        raise RuntimeError(f"{target.name}: expected {expected_bytes} bytes, received {total}")
    temp.replace(target)
    print(f"downloaded: {target.name} ({total:,} bytes)")


def main() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open(newline="") as f:
        for row in csv.DictReader(f):
            download(row["source_url"], INPUTS / row["file"], int(row["expected_bytes"]))
    download(CHAIN_URL, INPUTS / "hg18ToHg19.over.chain.gz")


if __name__ == "__main__":
    main()
