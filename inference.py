"""Run pretrained and retrained DeepSEA inference and collect Figure 2 arrays.

This orchestration script calls the project's two existing inference runners.
It does not train either model. Required model/data files and runner scripts are
located under the fuxian3 project root (see README.md).
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_PROJECT_ROOT = next(
    (candidate for candidate in (HERE, *HERE.parents)
     if (candidate / "run_inference.py").is_file()
     and (candidate / "run_inference_ours.py").is_file()),
    HERE,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT,
                        help="fuxian3 directory containing both run_inference scripts and model/data files")
    parser.add_argument("--device", default="cuda", choices=("cuda", "cpu"))
    parser.add_argument("--batch-size", type=int, default=512)
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()

    runners = {
        "pretrained": root / "run_inference.py",
        "ours": root / "run_inference_ours.py",
    }
    required = [
        *runners.values(),
        root / "deepsea_train" / "test.mat",
        root / "newdata" / "41592_2015_BFnmeth3547_MOESM647_ESM.csv",
        root / "reference_genome",
        root / "models" / "deepsea_predict.pth",
        root / "models" / "deepsea_variant_effects.pth",
        root / "training_checkpoints" / "best_model_FINAL_EPOCH53.pth",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        parser.error("missing inference inputs:\n  " + "\n  ".join(missing))
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")

    for name, runner in runners.items():
        print(f"\n=== {name}: {runner.name} ===", flush=True)
        subprocess.run(
            [sys.executable, str(runner), "--all", "--device", args.device,
             "--batch_size", str(args.batch_size)],
            cwd=root,
            check=True,
        )


    outputs = {
        "test_labels.npy": root / "test_labels.npy",
        "test_predictions_pretrained.npy": root / "test_predictions.npy",
        "dgf_ref_preds_pretrained.npy": root / "reproduction_results/gpu/dgf_ref_preds.npy",
        "dgf_alt_preds_pretrained.npy": root / "reproduction_results/gpu/dgf_alt_preds.npy",
        "test_predictions_ours.npy": root / "results/ours/test_predictions_ours.npy",
        "dgf_ref_preds_ours.npy": root / "reproduction_results/gpu/dgf_ref_preds_ours.npy",
        "dgf_alt_preds_ours.npy": root / "reproduction_results/gpu/dgf_alt_preds_ours.npy",
    }
    absent = [str(src) for src in outputs.values() if not src.is_file()]
    if absent:
        raise FileNotFoundError("inference runners did not produce expected files:\n  " + "\n  ".join(absent))
    output_dir = HERE / "experiment" / "Figure2" / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, source in outputs.items():
        destination = output_dir / filename
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)
        print(f"Saved {destination} ({destination.stat().st_size:,} bytes)")

    print("Inference arrays collected. Run plot.py pretrained and plot.py ours to draw the two plots.")


if __name__ == "__main__":
    main()
