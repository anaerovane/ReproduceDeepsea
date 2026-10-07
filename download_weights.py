"""Download DeepSEA checkpoints published in the Hugging Face repository."""

import argparse
from pathlib import Path

from huggingface_hub import hf_hub_download

REPO_ID = "aer0vane/reproduce_deepsea"
WEIGHTS = {
    "ours": "training_checkpoints/best_model_FINAL_EPOCH53.pth",
    "pretrained": "models/deepsea_predict.pth",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--which", choices=("ours", "pretrained", "both"), default="both",
        help="which checkpoint to download (default: both)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent,
        help="destination root; Hugging Face subdirectories are preserved",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    selected = WEIGHTS if args.which == "both" else {args.which: WEIGHTS[args.which]}
    for label, filename in selected.items():
        path = hf_hub_download(
            repo_id=REPO_ID,
            filename=filename,
            local_dir=args.output_dir,
        )
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
