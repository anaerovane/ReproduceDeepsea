#!/usr/bin/env python3
"""Merge CHT outputs and apply within-mark Benjamini-Hochberg FDR."""
from pathlib import Path
import sys
import pandas as pd
import numpy as np


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1.0)
    return q

here = Path(__file__).resolve().parent
parts = []
for mark in ("H3K4me3", "H3K27ac"):
    p = here / "results" / f"{mark}_cht_results.tsv"
    d = pd.read_csv(p, sep="\t")
    d["MARK"] = mark
    d["P.VALUE"] = pd.to_numeric(d["P.VALUE"], errors="coerce")
    d = d[d["P.VALUE"].notna() & d["P.VALUE"].between(0, 1)].copy()
    d["FDR"] = bh_fdr(d["P.VALUE"].to_numpy())
    parts.append(d)
out = pd.concat(parts, ignore_index=True)
out.to_csv(here / "results" / "mcvicker_cht_qtls.tsv", sep="\t", index=False)
for mark, d in out.groupby("MARK"):
    sig = d[d.FDR < 0.1]
    print(f"{mark}: tested pairs={len(d):,}, FDR<0.1 pairs={len(sig):,}, unique SNPs={sig['TEST.SNP.POS'].nunique():,}")
