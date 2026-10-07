#!/usr/bin/env python3
"""Run the configured pretrained and ours checkpoints on Experiment 6 inputs."""
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CONFIG = json.loads((HERE / "experiment6_config.json").read_text())
BASES = CONFIG["base_order"]
BASE_INDEX = {base: index for index, base in enumerate(BASES)}


class DeepSEA_Ours(nn.Module):
    """Architecture matching the local epoch-53 checkpoint."""

    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(4, 320, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2),
            nn.Conv2d(320, 480, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2),
            nn.Conv2d(480, 960, (1, 8)), nn.ReLU(), nn.Dropout(0.5),
        )
        self.classifier = nn.Sequential(
            nn.Linear(50880, 925), nn.ReLU(),
            nn.Linear(925, CONFIG["num_output_features"]), nn.Sigmoid()
        )

    def forward(self, x):
        x = self.features(x)
        return self.classifier(x.reshape(x.size(0), -1))


class LambdaBase(nn.Sequential):
    def __init__(self, fn, *args):
        super().__init__(*args)
        self.lambda_func = fn

    def forward_prepare(self, value):
        outputs = [module(value) for module in self._modules.values()]
        return outputs if outputs else value


class Lambda(LambdaBase):
    def forward(self, value):
        return self.lambda_func(self.forward_prepare(value))


class ReCodeAlphabet(nn.Module):
    def forward(self, value):
        return torch.stack([value[:, index, ...] for index in [0, 2, 1, 3]], dim=1)


class ConcatenateRC(nn.Module):
    def forward(self, value):
        reverse_complement = value.flip(dims=[1, 3])
        return torch.cat([value, reverse_complement], dim=0)


class AverageRC(nn.Module):
    def forward(self, value):
        count = value.shape[0] // 2
        return (value[:count] + value[count:]) / 2


def build_pretrained_model():
    backbone = nn.Sequential(
        nn.Conv2d(4, 320, (1, 8)), nn.Threshold(0, 1e-6),
        nn.MaxPool2d((1, 4), (1, 4)), nn.Dropout(0.2),
        nn.Conv2d(320, 480, (1, 8)), nn.Threshold(0, 1e-6),
        nn.MaxPool2d((1, 4), (1, 4)), nn.Dropout(0.2),
        nn.Conv2d(480, 960, (1, 8)), nn.Threshold(0, 1e-6), nn.Dropout(0.5),
        Lambda(lambda x: x.reshape(x.size(0), -1)),
        nn.Sequential(Lambda(lambda x: x.reshape(1, -1) if x.ndim == 1 else x),
                      nn.Linear(50880, 925)),
        nn.Threshold(0, 1e-6),
        nn.Sequential(Lambda(lambda x: x.reshape(1, -1) if x.ndim == 1 else x),
                      nn.Linear(925, CONFIG["num_output_features"])),
        nn.Sigmoid(),
    )
    return nn.Sequential(ReCodeAlphabet(), ConcatenateRC(), backbone, AverageRC())


def decode_onehot(one_hot):
    if one_hot.ndim != 2 or one_hot.shape[0] != len(BASES):
        raise ValueError(f"Expected ({len(BASES)}, sequence_length) one-hot input")
    if not np.allclose(one_hot.sum(axis=0), 1.0) or not np.isfinite(one_hot).all():
        raise ValueError("Sequence input must contain one unambiguous base at every position")
    return ''.join(BASES[index] for index in one_hot.argmax(axis=0))


def generate_mutations(sequence):
    original = decode_onehot(sequence)
    mutations, mutation_info = [], []
    for position, reference_base in enumerate(original):
        for alternate_base in BASES:
            if alternate_base == reference_base:
                continue
            mutated = sequence.copy()
            mutated[:, position] = 0
            mutated[BASE_INDEX[alternate_base], position] = 1.0
            mutations.append(mutated)
            mutation_info.append((position, reference_base, alternate_base))
    return np.asarray(mutations), mutation_info


def log2_odds_change(reference_probabilities, mutant_probabilities):
    epsilon = CONFIG["odds_probability_clip_epsilon"]
    reference = np.clip(reference_probabilities, epsilon, 1 - epsilon)
    mutant = np.clip(mutant_probabilities, epsilon, 1 - epsilon)
    reference_odds = reference / (1 - reference)
    mutant_odds = mutant / (1 - mutant)
    return np.log2(mutant_odds / reference_odds)


def load_model(model_config, device):
    if model_config["id"] == "ours":
        model = DeepSEA_Ours()
        checkpoint = torch.load(PROJECT / model_config["checkpoint"],
                                map_location=device, weights_only=False)
        state_dict = checkpoint.get("model_state_dict", checkpoint)
    elif model_config["id"] == "pretrained":
        model = build_pretrained_model()
        state_dict = torch.load(PROJECT / model_config["checkpoint"], map_location=device)
    else:
        raise ValueError(f"Unsupported model id: {model_config['id']}")
    model.load_state_dict(state_dict, strict=True)
    return model.to(device).eval()


def predict(model, sequences, model_config, device):
    channel_order = model_config["input_channel_order"]
    inputs = sequences[:, channel_order, :]
    batch_size = CONFIG["inference_batch_size"]
    probabilities = []
    with torch.no_grad():
        for start in range(0, len(inputs), batch_size):
            batch = torch.as_tensor(inputs[start:start + batch_size], dtype=torch.float32)
            batch = batch.unsqueeze(2).to(device)
            if model_config["reverse_complement_averaging"] == "inference_script":
                prediction = (model(batch) + model(batch.flip(dims=[1, 3]))) / 2
            elif model_config["reverse_complement_averaging"] == "checkpoint_internal":
                prediction = model(batch)
            else:
                raise ValueError("Unknown reverse-complement averaging setting")
            probabilities.append(prediction.cpu().numpy())
    return np.concatenate(probabilities, axis=0)


def main():
    torch.set_num_threads(4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)

    variants = CONFIG["variants"]
    sequences = {}
    for variant in variants:
        path = HERE / variant["sequence"]
        sequence = np.load(path, allow_pickle=False)
        expected_shape = (len(BASES), CONFIG["sequence_length"])
        if sequence.shape != expected_shape:
            raise ValueError(f"{path}: expected {expected_shape}, got {sequence.shape}")
        decoded = decode_onehot(sequence)
        if decoded[CONFIG["center_index"]] != variant["ref"]:
            raise ValueError(f"Reference base mismatch at {variant['chrom']}:{variant['position']}")
        sequences[(variant["chrom"], variant["position"])] = sequence

    effects, raw_predictions = {}, {}
    for model_config in CONFIG["models"]:
        model_id = model_config["id"]
        model = load_model(model_config, device)
        print(f"Model: {model_config['display_name']}", flush=True)
        for variant in variants:
            chrom, position = variant["chrom"], variant["position"]
            key = f"{model_id}_{chrom}_{position}"
            sequence = sequences[(chrom, position)]
            mutations, _ = generate_mutations(sequence)
            if len(mutations) != (len(BASES) - 1) * CONFIG["sequence_length"]:
                raise ValueError(f"Unexpected mutation count at {chrom}:{position}")

            original_probability = predict(model, sequence[None, ...], model_config, device)[0]
            mutant_probabilities = predict(model, mutations, model_config, device)
            effect = log2_odds_change(original_probability, mutant_probabilities)
            effects[key] = effect
            raw_predictions[f"{key}_original"] = original_probability
            raw_predictions[f"{key}_mutated"] = mutant_probabilities
            print(f"  {chrom}:{position}: {len(mutations)} substitutions", flush=True)
        del model

    np.savez(HERE / "experiment6_all_data.npz", **effects)
    np.savez_compressed(HERE / "experiment6_raw_predictions.npz", **raw_predictions)
    print("Inference complete; generating figures from the saved effects.", flush=True)
    from experiment6_plot import main as plot_main
    plot_main()


if __name__ == "__main__":
    main()
