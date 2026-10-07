#!/usr/bin/env python3
"""Plot configured Experiment 6 effects from the saved inference cache."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CONFIG = json.loads((HERE / "experiment6_config.json").read_text())
BASES = CONFIG["base_order"]
CMAP = LinearSegmentedColormap.from_list(
    "blue_white_yellow", ["#0571b0", "#ffffff", "#ffff00"], N=256
)


def effect_matrix(values, sequence, feature_index):
    sequence_length = CONFIG["sequence_length"]
    expected_shape = (
        (len(BASES) - 1) * sequence_length,
        CONFIG["num_output_features"],
    )
    if values.shape != expected_shape or not np.isfinite(values).all():
        raise ValueError(f"Expected finite {expected_shape} effects")
    if sequence.shape != (len(BASES), sequence_length):
        raise ValueError("Unexpected one-hot sequence shape")
    if not np.allclose(sequence.sum(axis=0), 1.0):
        raise ValueError("Sequence has ambiguous or malformed one-hot columns")

    original = ''.join(BASES[index] for index in sequence.argmax(axis=0))
    effects = np.zeros((len(BASES), sequence_length), dtype=np.float64)
    mutation_index = 0
    for position, reference_base in enumerate(original):
        for row, alternate_base in enumerate(BASES):
            if alternate_base == reference_base:
                continue
            effects[row, position] = values[mutation_index, feature_index]
            mutation_index += 1
    if mutation_index != len(values):
        raise ValueError("Mutation order does not match the configured sequence")
    return effects, original


def main():
    names_path = PROJECT / CONFIG["predictor_names"]
    predictor_names = names_path.read_text().splitlines()
    with np.load(HERE / "experiment6_all_data.npz", allow_pickle=False) as data:
        panels = []
        metadata = {
            "mutation_order": CONFIG["mutation_order"],
            "effect": CONFIG["effect"],
            "panels": [],
        }
        for variant in CONFIG["variants"]:
            chrom, position = variant["chrom"], variant["position"]
            feature_index = variant["feature_index"]
            if feature_index >= len(predictor_names):
                raise ValueError(f"Feature index is out of range: {feature_index}")
            actual_name = predictor_names[feature_index]
            if actual_name != variant["feature_name"]:
                raise ValueError(
                    f"Feature mismatch at index {feature_index}: "
                    f"expected {variant['feature_name']!r}, got {actual_name!r}"
                )

            sequence = np.load(HERE / variant["sequence"], allow_pickle=False)
            reference = ''.join(BASES[index] for index in sequence.argmax(axis=0))
            center = CONFIG["center_index"]
            if reference[center] != variant["ref"]:
                raise ValueError(f"Reference mismatch at {chrom}:{position}")

            matrices = {}
            for model in CONFIG["models"]:
                key = f"{model['id']}_{chrom}_{position}"
                matrices[model["id"]], _ = effect_matrix(
                    data[key], sequence, feature_index
                )
            color_limit = max(float(np.abs(matrix).max()) for matrix in matrices.values())
            color_limit = max(color_limit, 1e-6)
            panels.append((variant, matrices, reference, color_limit))
            metadata["panels"].append({
                "tf": variant["tf"],
                "cell": variant["cell"],
                "chrom": chrom,
                "position": position,
                "feature_index": feature_index,
                "feature_name": actual_name,
                "color_limit": color_limit,
            })

        for model in CONFIG["models"]:
            model_id = model["id"]
            figure = plt.figure(figsize=(9, 15), facecolor="white")
            grid = figure.add_gridspec(
                len(panels), 1, left=.045, right=.96, bottom=.035, top=.96, hspace=.27
            )
            figure.suptitle(model["display_name"], fontsize=14, y=.99)
            for row, (variant, matrices, reference, color_limit) in enumerate(panels):
                subgrid = grid[row].subgridspec(
                    2, 2, width_ratios=[100, 1], hspace=.24, wspace=.055
                )
                zoom_start, zoom_end = CONFIG["zoom_interval"]
                for column, (start, end) in enumerate(
                    ((0, CONFIG["sequence_length"]), (zoom_start, zoom_end))
                ):
                    axis = figure.add_subplot(subgrid[column, 0])
                    color_axis = figure.add_subplot(subgrid[column, 1])
                    image = axis.imshow(
                        matrices[model_id][:, start:end], origin="upper",
                        extent=(start, end, len(BASES), 0), aspect="auto",
                        cmap=CMAP, vmin=-color_limit, vmax=color_limit,
                        interpolation="none",
                    )
                    axis.set_xlim(start, end)
                    axis.set_yticks([])
                    axis.set_xticks(
                        np.arange(0, CONFIG["sequence_length"] + 1, 100)
                        if column == 0 else [zoom_start, CONFIG["center_index"], zoom_end]
                    )
                    axis.xaxis.tick_top()
                    axis.tick_params(axis="x", labelsize=5, length=2, pad=1, color="#888888")
                    for spine in axis.spines.values():
                        spine.set_color("#888888")
                        spine.set_linewidth(.6)
                    if column == 0:
                        axis.set_title(
                            f"{variant['tf']} ({variant['cell']})",
                            loc="left", fontsize=12, pad=12,
                        )
                    else:
                        for index in range(start, end):
                            axis.text(
                                index + .5, -.025, reference[index],
                                transform=axis.get_xaxis_transform(), ha="center",
                                va="top", fontsize=3.1, family="monospace", clip_on=False,
                            )
                    colorbar = figure.colorbar(image, cax=color_axis)
                    colorbar.set_ticks(np.linspace(-color_limit, color_limit, 9))
                    colorbar.ax.tick_params(labelsize=4, length=1, pad=1)
                    colorbar.outline.set_linewidth(.4)

            for extension in ("png", "pdf"):
                output = HERE / f"experiment6_saturation_{model_id}.{extension}"
                figure.savefig(output, dpi=240)
                print(f"Saved: {output}", flush=True)
            plt.close(figure)

    (HERE / "experiment6_plot_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
