"""Fit the unchanged WASP AS-dispersion likelihood per donor in parallel."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import subprocess
import sys

here = Path(__file__).resolve().parent

wasp_dir = os.environ.get("WASP_CHT_DIR")
if not wasp_dir:
    raise SystemExit("Set WASP_CHT_DIR to the directory containing fit_as_coefficients_vectorized.py")
wasp = Path(wasp_dir) / "fit_as_coefficients_vectorized.py"
results = here / "results"
jobs = []
for mark in ("H3K4me3", "H3K27ac"):
    inputs = [x.strip() for x in (results / f"{mark}_input_files.txt").read_text().splitlines() if x.strip()]
    for i, file in enumerate(inputs):
        list_path = results / f"{mark}_as_one_{i:02d}.txt"
        out_path = results / f"{mark}_as_one_{i:02d}.dispersion"
        list_path.write_text(file + "\n")
        jobs.append((mark, i, list_path, out_path))

def run(job):
    mark, i, list_path, out_path = job
    p = subprocess.run([sys.executable, str(wasp), str(list_path), str(out_path)],
                       capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"{mark} donor {i} failed:\n{p.stderr[-4000:]}")
    return mark, i, p.stderr

with ThreadPoolExecutor(max_workers=16) as pool:
    futures = [pool.submit(run, job) for job in jobs]
    for f in as_completed(futures):
        mark, i, log = f.result()
        lines = [x for x in log.splitlines() if "dispersion estimate:" in x]
        print(f"{mark} donor {i}: {lines[-1] if lines else 'finished'}", flush=True)

for mark in ("H3K4me3", "H3K27ac"):
    values = [(results / f"{mark}_as_one_{i:02d}.dispersion").read_text().strip()
              for i in range(10)]
    (results / f"{mark}_as_dispersion.txt").write_text("\n".join(values) + "\n")
