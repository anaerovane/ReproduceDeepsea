"""Run allele-specific DeepSEA predictions and draw Supplementary Figure 4.

Input is the CHT result table produced from the McVicker GSE47991 read-count
files. A single row per SNP/mark is retained (the most significant target
region); the CHT alpha/beta estimates provide the observed REF/ALT direction.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import gzip
from bisect import bisect_left, bisect_right
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import run_inference as pretrained_inference
import run_inference_ours as ours_inference
from model_assets import (
    resolve_ours_checkpoint,
    resolve_pretrained_variant_effects_checkpoint,
)

MODEL_INFO = {
    "pretrained": {
        "path": resolve_pretrained_variant_effects_checkpoint(ROOT),
        "sha256": "98b20183d3288f220153a93d75bca1e6fb76800a931bf4f8b546382fdd8c99ea",
    },
    "ours": {
        "path": resolve_ours_checkpoint(ROOT),
        "sha256": "5979e5b48ef892b76f23aa3af4149787e8f7a0ef5f5bfef5e676997a6c1731bb",
    },
}
PREDICTORS = {"H3K4me3": 860, "H3K27ac": 855}
MARKS = ["H3K4me3", "H3K27ac"]
OUT = HERE / "results"
REF_DIR = ROOT / "reference_genome"
CHAIN = HERE / "inputs/hg18ToHg19.over.chain.gz"


def load_qtls(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"MARK", "TEST.SNP.CHROM", "TEST.SNP.POS", "TEST.SNP.ID",
                "TEST.SNP.REF.ALLELE", "TEST.SNP.ALT.ALLELE", "P.VALUE",
                "FDR", "ALPHA", "BETA"}
    miss = required - set(df.columns)
    if miss:
        raise ValueError(f"QTL file missing columns: {sorted(miss)}")
    df = df[df["MARK"].isin(MARKS)].copy()
    for col in ("TEST.SNP.POS", "P.VALUE", "FDR", "ALPHA", "BETA"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df["FDR"] < 0.1) & df["TEST.SNP.POS"].notna() &
            df["ALPHA"].notna() & df["BETA"].notna()].copy()


    df = df.sort_values(["MARK", "FDR", "P.VALUE"]).drop_duplicates(
        ["MARK", "TEST.SNP.CHROM", "TEST.SNP.POS", "TEST.SNP.REF.ALLELE", "TEST.SNP.ALT.ALLELE"])
    if df.empty:
        raise ValueError("No unique SNPs passed FDR < 0.1")
    df["truth_alt_high"] = df["BETA"] > df["ALPHA"]
    return df.reset_index(drop=True)


def one_hot(seq: str) -> np.ndarray:
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    out = np.zeros((4, len(seq)), dtype=np.float32)
    for i, b in enumerate(seq):
        out[idx[b], i] = 1.0
    return out


def lift_hg18_positions(df: pd.DataFrame, chain_path: Path):
    """Map hg18 source coordinates on chain t to hg19 destination coordinates on q."""
    wanted = {}
    for i, row in df.iterrows():
        chrom = str(row["TEST.SNP.CHROM"])
        wanted.setdefault(chrom, []).append((int(row["TEST.SNP.POS"]) - 1, i))
    for chrom in wanted:
        wanted[chrom].sort()
    positions = {c: [x[0] for x in entries] for c, entries in wanted.items()}
    mapped = {}
    with gzip.open(chain_path, "rt") as fh:
        for line in fh:
            if not line.startswith("chain "):
                continue
            h = line.split()
            score = int(h[1])


            tname, tpos = h[2], int(h[5])
            qname, qsize, qstrand, qpos = h[7], int(h[8]), h[9], int(h[10])
            while True:
                fields = fh.readline().split()
                if not fields:
                    break
                size = int(fields[0])
                thi = tpos + size
                if tname in positions:
                    arr = positions[tname]
                    lo, hi = bisect_left(arr, tpos), bisect_left(arr, thi)
                    for tbase, idx in wanted[tname][lo:hi]:
                        delta = tbase - tpos
                        target0 = qpos + delta if qstrand == "+" else qsize - 1 - (qpos + delta)
                        old = mapped.get(idx)
                        if old is None or score > old[0]:
                            mapped[idx] = (score, qname, target0, qstrand)
                if len(fields) == 1:
                    break
                tpos += size + int(fields[1])
                qpos += size + int(fields[2])
    return mapped


def fetch_variant_sequences(df: pd.DataFrame, chain_path: Path = CHAIN):
    genomes = {}
    refs, alts, keep, orientation = [], [], [], []
    complement = str.maketrans("ACGT", "TGCA")
    lifted = lift_hg18_positions(df, chain_path)
    for ix, row in df.iterrows():
        mapping = lifted.get(ix)
        if mapping is None:
            continue
        _, chrom, pos, chain_strand = mapping
        if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
            continue
        if chrom not in genomes:
            p = REF_DIR / f"{chrom}.fa"
            genomes[chrom] = "".join(x.strip() for x in p.open() if not x.startswith(">" )).upper() if p.exists() else None
        genome = genomes[chrom]
        if genome is None:
            continue
        ref, alt = str(row["TEST.SNP.REF.ALLELE"]).upper(), str(row["TEST.SNP.ALT.ALLELE"]).upper()


        if chain_strand == "-":
            ref, alt = ref.translate(complement), alt.translate(complement)
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            continue


        start, end = pos - 499, pos + 501
        if start < 0 or end > len(genome):
            continue
        seq = genome[start:end]
        if len(seq) != 1000 or any(b not in "ACGT" for b in seq):
            continue
        genome_ref = seq[499]
        direct_ref, direct_alt = genome_ref == ref, genome_ref == alt
        comp_ref = genome_ref == ref.translate(complement)
        comp_alt = genome_ref == alt.translate(complement)
        matches = [direct_ref, direct_alt, comp_ref, comp_alt]
        if sum(matches) != 1:


            continue
        if direct_ref:
            genome_alt, flip_truth = alt, False
            orientation_label = "direct_ref_alt"
        elif direct_alt:
            genome_alt, flip_truth = ref, True
            orientation_label = "direct_swapped"
        elif comp_ref:
            genome_alt, flip_truth = alt.translate(complement), False
            orientation_label = "complement_ref_alt"
        else:
            genome_alt, flip_truth = ref.translate(complement), True
            orientation_label = "complement_swapped"
        alt_seq = seq[:499] + genome_alt + seq[500:]
        refs.append(one_hot(seq)); alts.append(one_hot(alt_seq)); keep.append(ix)
        orientation.append((orientation_label, flip_truth, genome_ref, genome_alt, chrom, pos + 1,
                            chain_strand, int(mapping[0])))
    return np.asarray(refs), np.asarray(alts), np.asarray(keep, dtype=int), orientation


def run_predictions(refs: np.ndarray, alts: np.ndarray, device: str, batch_size: int):
    pre = pretrained_inference.load_variant_effect_model(device)
    ours = ours_inference.load_our_model(device)
    outputs = {}
    for label, model in [("pretrained", pre), ("ours", ours)]:
        pred_ref, pred_alt = [], []


        if label == "ours":
            model_refs = refs[:, [0, 2, 1, 3], :]
            model_alts = alts[:, [0, 2, 1, 3], :]
        else:
            model_refs, model_alts = refs, alts
        with torch.inference_mode():
            for i in range(0, len(model_refs), batch_size):
                xr = torch.from_numpy(model_refs[i:i + batch_size]).unsqueeze(2).to(device)
                xa = torch.from_numpy(model_alts[i:i + batch_size]).unsqueeze(2).to(device)
                if label == "pretrained":


                    pr = (model(xr) + model(xr.flip(dims=[1, 3]))) / 2
                    pa = (model(xa) + model(xa.flip(dims=[1, 3]))) / 2
                else:
                    pr = (model(xr) + model(xr.flip(dims=[1, 3]))) / 2
                    pa = (model(xa) + model(xa.flip(dims=[1, 3]))) / 2
                pred_ref.append(pr.cpu().numpy()); pred_alt.append(pa.cpu().numpy())
        outputs[label] = (np.concatenate(pred_ref), np.concatenate(pred_alt))
        del model
    return outputs


def write_figures(df: pd.DataFrame, raw: dict, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    curve_rows = []
    for model_name, (p_ref, p_alt) in raw.items():
        d = df.copy()
        d["p_ref"] = [p_ref[i, PREDICTORS[m]] for i, m in enumerate(d["MARK"])]
        d["p_alt"] = [p_alt[i, PREDICTORS[m]] for i, m in enumerate(d["MARK"])]
        d["signed_delta"] = d["p_alt"] - d["p_ref"]
        d["margin"] = d["signed_delta"].abs()
        d["pred_alt_high"] = d["signed_delta"] > 0
        d["correct"] = d["pred_alt_high"] == d["truth_alt_high"]
        d.to_csv(output_dir / f"{model_name}_variant_predictions.tsv", sep="\t", index=False)

        fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), sharey=True)
        for ax, mark, xmax in zip(axes, MARKS, [0.06, 0.08]):
            x = np.linspace(0, xmax, 301)
            sub = d[d["MARK"] == mark]
            y, n = [], []
            for threshold in x:
                selected = sub[sub["margin"] >= threshold]
                n.append(len(selected))
                accuracy = float(selected["correct"].mean()) if len(selected) else np.nan
                y.append(accuracy)
                curve_rows.append({"model": model_name, "mark": mark, "margin_threshold": threshold,
                                   "accuracy": accuracy, "n_above_threshold": len(selected),
                                   "n_qtl_snps": len(sub)})
            ax.step(x, y, where="post", color="#3977ad", linewidth=1.2)
            ax.set_xlim(0, xmax); ax.set_ylim(0.3, 1.0)
            ax.set_xlabel("Margin")
            ax.set_title(mark)
            ax.grid(False)
            ax.text(0.98, 0.04, f"n={len(sub)} QTL SNPs", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=8, color="#555555")
        axes[0].set_ylabel("Accuracy above margin threshold")
        fig.suptitle(f"Diagnostic reconstruction — {model_name}", y=1.02)
        fig.tight_layout()
        fig.savefig(output_dir / f"supplementary_figure4_{model_name}.png", dpi=300, bbox_inches="tight")
        fig.savefig(output_dir / f"supplementary_figure4_{model_name}.pdf", bbox_inches="tight")
        plt.close(fig)

    pd.DataFrame(curve_rows).to_csv(output_dir / "accuracy_vs_margin.tsv", sep="\t", index=False)

    np.savez_compressed(output_dir / "histone_qtl_raw_predictions.npz",
                        pred_ref_pretrained=raw["pretrained"][0], pred_alt_pretrained=raw["pretrained"][1],
                        pred_ref_ours=raw["ours"][0], pred_alt_ours=raw["ours"][1],
                        qtl_row_index=df.index.to_numpy())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qtl-table", type=Path, default=HERE / "results/mcvicker_cht_qtls.tsv")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--chain", type=Path, default=CHAIN)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    df = load_qtls(args.qtl_table)
    qtl_candidate_count = len(df)
    refs, alts, keep, orientation = fetch_variant_sequences(df, args.chain)
    df = df.iloc[keep].reset_index(drop=True)
    df["allele_orientation"] = [x[0] for x in orientation]
    df["genome_ref"] = [x[2] for x in orientation]
    df["genome_alt"] = [x[3] for x in orientation]
    df["hg19_chrom"] = [x[4] for x in orientation]
    df["hg19_pos"] = [x[5] for x in orientation]
    df["chain_strand"] = [x[6] for x in orientation]
    df["lift_chain_score"] = [x[7] for x in orientation]
    flips = np.asarray([x[1] for x in orientation], dtype=bool)
    df["truth_alt_high"] = np.logical_xor(df["truth_alt_high"].to_numpy(dtype=bool), flips)
    print(f"QTL SNPs passing FDR < 0.1 and allele harmonization to hg19: {len(df)}/{qtl_candidate_count}")
    if not len(df):
        raise SystemExit("No usable variant sequences")
    predictions = run_predictions(refs, alts, args.device, args.batch_size)
    write_figures(df, predictions, args.output_dir)
    counts = {str(k): int(v) for k, v in df.groupby("MARK").size().to_dict().items()}
    report = {
        "cohort_source": "McVicker et al. 2013, GEO GSE47991 processed read counts",
        "source_build": "hg18 (GEO GSE47991 processed read counts)",
        "reference_build": "hg19 / GRCh37",
        "variant_center_index_0based": 499,
        "liftover": {"chain_file": str(args.chain.relative_to(HERE)), "direction": "hg18 to hg19", "method": "map source hg18 positions on UCSC chain t side to destination hg19 q side; highest-score primary chromosome mapping retained"},
        "qtl_method": "WASP combined haplotype test; FDR < 0.1; one most significant region per unique SNP and mark",
        "wasp_source_commit": "d3b8447fd7719fffa00b856fd1f27c845554693e",
        "fdr_method": "Benjamini-Hochberg, independently per histone mark",
        "deepsea_outputs": {"H3K4me3": 860, "H3K27ac": 855},
        "pretrained_inference_mode": "deepsea_variant_effects.pth with forward/reverse-complement mean, matching the DeepSEA variantEffects workflow",
        "input_channel_order": {
            "pretrained": "ACGT one-hot; variant-effect wrapper remaps to AGCT; forward and reverse-complement predictions averaged",
            "ours": "AGCT one-hot",
        },
        "allele_harmonization": {
            "rule": "map hg18 chain t coordinates to hg19 q coordinates; orient alleles for negative chain blocks; match CHT allele pairs to hg19 reference and swap effect direction if ALT is the reference base; exclude mismatches and strand-ambiguous palindromic SNPs",
            "orientation_counts": {str(k): int(v) for k, v in df["allele_orientation"].value_counts().to_dict().items()},
            "unharmonized_rows_excluded": int(qtl_candidate_count - len(df)),
        },
        "models": {k: {"path": str(v["path"].relative_to(ROOT)), "sha256": v["sha256"]}
                   for k, v in MODEL_INFO.items()},
        "qtl_counts_after_allele_harmonization": counts,
        "predicted_count": len(df),
        "device": args.device,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
