#!/usr/bin/env python3
"""Re-run experiment 7 ref/alt predictions from local FASTA and checkpoints."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint

OUT = Path(__file__).resolve().parent
FASTA = ROOT / "reference_genome"
PREDICTORS = (ROOT / "predictor_names.txt").read_text().splitlines()
CHECKPOINTS = {
    "ours_epoch53": resolve_ours_checkpoint(ROOT),
    "pretrained": resolve_pretrained_predict_checkpoint(ROOT),
}
VARIANTS = [
    {
        "id": "rs4784227",
        "chrom": "chr16",
        "pos_1based": 52599188,
        "ref": "C",
        "alt": "T",
        "reported_factor": "FOXA1",
        "local_manifest_factor": "FOXA2",
    },
    {
        "id": "alpha_thalassemia_chr16_209709",
        "chrom": "chr16",
        "pos_1based": 209709,
        "ref": "T",
        "alt": "C",
        "reported_factor": "GATA1",
        "local_manifest_factor": "GATA1",
    },
    {
        "id": "pancreatic_agenesis_chr10_23508363",
        "chrom": "chr10",
        "pos_1based": 23508363,
        "ref": "A",
        "alt": "G",
        "reported_factor": "FOXA2",
        "local_manifest_factor": "FOXA1",
    },
]


def load_architecture():
    source = ROOT / "experiment/Supplementary_Figure3/experiment6_saturation_full.py"
    spec = importlib.util.spec_from_file_location("experiment6_architecture", source)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def load_sequence(variant: dict) -> str:
    path = FASTA / f"{variant['chrom']}.fa"
    sequence = "".join(line.strip() for line in path.open() if not line.startswith(">"))
    sequence = sequence.upper()
    center = variant["pos_1based"] - 1
    window = sequence[center - 500:center + 500]
    if len(window) != 1000:
        raise ValueError(f"Expected 1000 bases for {variant['id']}; got {len(window)}")
    observed_ref = window[500]
    if observed_ref != variant["ref"]:
        raise ValueError(
            f"Reference mismatch at {variant['chrom']}:{variant['pos_1based']}: "
            f"expected {variant['ref']}, FASTA has {observed_ref}"
        )
    if set(window) - set("AGCT"):
        raise ValueError(f"Ambiguous bases in {variant['id']} reference window")
    return window


def one_hot(sequence: str) -> np.ndarray:
    # Match experiment6_saturation_full.py exactly: channel order AGCT.
    mapping = {base: index for index, base in enumerate("AGCT")}
    result = np.zeros((4, len(sequence)), dtype=np.float32)
    for i, base in enumerate(sequence):
        result[mapping[base], i] = 1.0
    return result


def target_indices(factor: str) -> list[int]:
    if factor == "GATA1":
        names = [(i, name) for i, name in enumerate(PREDICTORS) if "GATA-1" in name.upper()]
    else:
        names = [(i, name) for i, name in enumerate(PREDICTORS) if factor in name.upper()]
    if not names:
        raise ValueError(f"No predictor index found for {factor}")
    # Preserve all matching cell/assay outputs: the report's FOXA1 probability
    # for rs4784227 is the T-47D assay, not the HepG2 assay in the old plot.
    return [i for i, _ in names]


def predict(model, sequence: str, model_name: str, device: torch.device) -> np.ndarray:
    encoded = one_hot(sequence)
    if model_name == "pretrained":
        # The local pretrained inference path consumes ACGT channel order.
        encoded = encoded[[0, 2, 1, 3], :].copy()
    tensor = torch.from_numpy(encoded).unsqueeze(0).unsqueeze(2).to(device)
    with torch.inference_mode():
        if model_name == "ours_epoch53":
            pred = (model(tensor) + model(tensor.flip(dims=[1, 3]))) / 2
        else:
            pred = model(tensor)
    return pred.detach().float().cpu().numpy()[0]


def main():
    torch.set_num_threads(4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    architecture = load_architecture()
    outputs = {}
    metadata = {
        "device": str(device),
        "reference": "local reference_genome FASTA; coordinates interpreted as 1-based",
        "window": "1000 bp, [pos-500, pos+499] around 1-based SNP position",
        "models": {},
        "variants": [],
        "outputs": {},
    }
    for name, path in CHECKPOINTS.items():
        metadata["models"][name] = {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        if name == "ours_epoch53":
            model = architecture.DeepSEA_Ours()
            checkpoint = torch.load(path, map_location=device, weights_only=False)
            state = checkpoint.get("model_state_dict", checkpoint) if isinstance(checkpoint, dict) else checkpoint
        else:
            model = architecture.build_predict_model()
            state = torch.load(path, map_location=device, weights_only=False)
        model.load_state_dict(state)
        model.to(device).eval()

        for variant in VARIANTS:
            ref_sequence = load_sequence(variant)
            alt_sequence = ref_sequence[:500] + variant["alt"] + ref_sequence[501:]
            p_ref = predict(model, ref_sequence, name, device)
            p_alt = predict(model, alt_sequence, name, device)
            outputs[f"{name}_{variant['id']}_ref"] = p_ref
            outputs[f"{name}_{variant['id']}_alt"] = p_alt
            factors = sorted({variant["reported_factor"], variant["local_manifest_factor"]})
            factor_values = {}
            for factor in factors:
                values = []
                for index in target_indices(factor):
                    values.append({
                        "index": index,
                        "predictor": PREDICTORS[index],
                        "ref_probability": float(p_ref[index]),
                        "alt_probability": float(p_alt[index]),
                        "delta": float(p_alt[index] - p_ref[index]),
                        "direction": "increase" if p_alt[index] > p_ref[index] else "decrease" if p_alt[index] < p_ref[index] else "unchanged",
                    })
                factor_values[factor] = values
            metadata["variants"].append({
                **variant,
                "observed_reference": ref_sequence[500],
                "window_ref": ref_sequence,
                "window_alt": alt_sequence,
                "factor_predictions": {name: factor_values},
            })
    # Consolidate per-model records under one entry per variant.
    by_variant = {}
    for record in metadata["variants"]:
        variant_id = record["id"]
        if variant_id not in by_variant:
            by_variant[variant_id] = {**record, "factor_predictions": {}}
        by_variant[variant_id]["factor_predictions"].update(record["factor_predictions"])
    metadata["variants"] = list(by_variant.values())
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "experiment7_raw_predictions.npz", **outputs)
    (OUT / "experiment7_audit.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    for variant in metadata["variants"]:
        print(variant["id"], variant["chrom"], variant["pos_1based"], variant["ref"] + ">" + variant["alt"])
        for model_name, prediction in variant["factor_predictions"].items():
            for factor, values in prediction.items():
                for item in values:
                    print(model_name, factor, item["predictor"],
                          f"ref={item['ref_probability']:.6f}",
                          f"alt={item['alt_probability']:.6f}", item["direction"])
    print(f"Saved raw predictions and audit to {OUT}")


if __name__ == "__main__":
    main()
