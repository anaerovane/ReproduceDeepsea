"""Extract and validate the four 1 kb GRCh37 reference windows for Experiment 6."""
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CONFIG = json.loads((HERE / "experiment6_config.json").read_text())
BASES = CONFIG["base_order"]
BASE_INDEX = {base: index for index, base in enumerate(BASES)}


def read_fasta(path):
    sequence = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                sequence.append(line.strip().upper())
    return "".join(sequence)


def main():
    length = CONFIG["sequence_length"]
    fasta_dir = PROJECT / CONFIG["sequence_source"]["fasta_dir"]
    for variant in CONFIG["variants"]:
        chrom = variant["chrom"]
        position = variant["position"]
        center = position - 1
        center_index = CONFIG["center_index"]
        start = center - center_index
        end = start + length
        genome = read_fasta(fasta_dir / f"{chrom}.fa")
        if start < 0 or end > len(genome):
            raise ValueError(f"Window is out of bounds: {chrom}:{position}")
        actual_ref = genome[center]
        if actual_ref != variant["ref"]:
            raise ValueError(
                f"Reference mismatch at {chrom}:{position}: "
                f"FASTA={actual_ref}, config={variant['ref']}"
            )
        sequence = genome[start:end]
        if len(sequence) != length or any(base not in BASE_INDEX for base in sequence):
            raise ValueError(f"Window has wrong length or ambiguous bases: {chrom}:{position}")

        one_hot = np.zeros((len(BASES), length), dtype=np.float32)
        for column, base in enumerate(sequence):
            one_hot[BASE_INDEX[base], column] = 1.0
        out = HERE / variant["sequence"]
        np.save(out, one_hot, allow_pickle=False)
        print(f"Saved {chrom}:{position} {variant['ref']}>{variant['alt']} -> {out}")


if __name__ == "__main__":
    main()
