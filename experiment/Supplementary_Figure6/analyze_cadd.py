"""Score all returned CADD annotations with empirical DeepSEA E-values.

Uses the existing pretrained and ours effect matrices from Figure3 and
cohort-specific random-negative SNPs as empirical backgrounds. Missing CADD
rows or fields remain missing and are reported in coverage and provenance.
"""
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"
FIGURE3 = HERE.parent / "Figure3" / "outputs"
N_CHROMATIN = 919
CADD_COLUMNS = ("priPhCons", "priPhyloP", "GerpN", "GerpS")
COHORTS = ("eqtl", "gwas")
NEGATIVE_GROUPS = (
    "Random negative SNP", "31kbp negative SNP", "6.3kbp negative SNP",
    "710bp negative SNP", "360bp negative SNP",
)


def canonical_chrom(s):
    return str(s).removeprefix("chr")


def read_cadd(path):
    with gzip.open(path, "rt") as handle:
        metadata_lines = 0
        for line in handle:
            if line.startswith("##"):
                metadata_lines += 1
            else:
                break
    usecols = ["#Chrom", "Pos", "Ref", "Alt", *CADD_COLUMNS]
    frame = pd.read_csv(path, sep="\t", compression="gzip", skiprows=metadata_lines,
                        usecols=usecols, low_memory=False)
    frame = frame.rename(columns={"#Chrom": "chrom", "Pos": "pos", "Ref": "ref", "Alt": "alt"})
    frame["chrom"] = frame.chrom.map(canonical_chrom)
    frame["pos"] = pd.to_numeric(frame.pos, errors="raise").astype(np.int64)
    frame["ref"] = frame.ref.str.upper()
    frame["alt"] = frame.alt.str.upper()
    for col in CADD_COLUMNS:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def load_annotations():
    files = sorted(OUTPUTS.glob("cadd_inclAnno_batch_*.tsv.gz"))
    if not files:
        raise RuntimeError("No downloaded CADD annotation batches found")
    raw = pd.concat((read_cadd(path) for path in files), ignore_index=True)
    keys = ["chrom", "pos", "ref", "alt"]
    duplicated = raw[raw.duplicated(keys, keep=False)]
    if not duplicated.empty:
        score_counts = duplicated.groupby(keys, dropna=False)[list(CADD_COLUMNS)].nunique(dropna=False)
        if (score_counts > 1).any(axis=None):
            raise RuntimeError("Overlapping CADD batches contain conflicting annotation values")
    raw = raw.drop_duplicates(keys, keep="first")

    variants = pd.read_csv(FIGURE3 / "variants.tsv", sep="\t",
                           usecols=["chr", "pos", "ref", "alt", "variant_id"])
    variants["chrom"] = variants.chr.map(canonical_chrom)
    variants["pos"] = variants.pos.astype(np.int64)
    variants["ref"] = variants.ref.str.upper()
    variants["alt"] = variants.alt.str.upper()
    variants = variants[["chrom", "pos", "ref", "alt", "variant_id"]]
    if variants.duplicated(["chrom", "pos", "ref", "alt"]).any():
        raise RuntimeError("Experiment8 variant table has duplicate variant keys")
    joined = raw.merge(variants, on=keys, how="left", validate="one_to_one", indicator=True)
    if not joined._merge.eq("both").all():
        raise RuntimeError(f"{int((joined._merge != 'both').sum())} CADD records did not map to Figure3 variants")
    return files, joined.drop(columns="_merge"), int(variants.variant_id.max()) + 1


def empirical_evalue(background, query):
    """Smoothed empirical upper-tail E-value, including exact ties as not greater."""
    bg = np.asarray(background, dtype=np.float64)
    q = np.asarray(query, dtype=np.float64)
    bg = bg[np.isfinite(bg)]
    result = np.full(q.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(q)
    if len(bg) == 0 or not valid.any():
        return result
    bg.sort()


    n_at_or_below = np.searchsorted(bg, q[valid], side="right")
    result[valid] = (len(bg) - n_at_or_below + 1) / (len(bg) + 1)
    return result


def load_cohorts(cadd_by_id, n_variants):
    cohorts = {}
    for name in COHORTS:
        path = FIGURE3 / f"{name}_pretrained_oof.tsv.gz"
        frame = pd.read_csv(path, sep="\t")
        ids = frame.variant_id.to_numpy(dtype=np.int64)
        if ids.min() < 0 or ids.max() >= n_variants:
            raise RuntimeError(f"{name}: variant_id outside effect matrix")
        frame["_vid"] = ids
        for idx, col in enumerate(CADD_COLUMNS):
            frame[col] = cadd_by_id[ids, idx]
        cohorts[name] = frame
    return cohorts


def prepare_effect_magnitude(path, n_variants):
    effects = np.load(path, mmap_mode="r")
    if effects.shape != (n_variants, 2 * N_CHROMATIN):
        raise RuntimeError(f"Unexpected effect matrix shape for {path.name}: {effects.shape}")


    magnitude = np.empty((n_variants, N_CHROMATIN), dtype=np.float32)
    for start in range(0, n_variants, 20_000):
        stop = min(start + 20_000, n_variants)
        block = effects[start:stop]
        magnitude[start:stop] = np.abs(block[:, :N_CHROMATIN]) * np.abs(block[:, N_CHROMATIN:])
        if start == 0 or stop == n_variants or (start // 20_000) % 20 == 0:
            print(f"{path.stem}: prepared effect magnitudes {stop:,}/{n_variants:,}", flush=True)
    return magnitude


def score_cohort(frame, magnitude, model, cohort, cached_chrom=None):
    ids = frame._vid.to_numpy(dtype=np.int64)
    bg_ids = np.unique(frame.loc[frame.label.eq("Random negative SNP"), "_vid"].to_numpy(dtype=np.int64))
    if len(bg_ids) == 0:
        raise RuntimeError(f"{cohort}: no random-negative background")
    if cached_chrom is not None:
        chrom_gm = cached_chrom
        print(f"{cohort}/{model}: reused cached chromatin scores", flush=True)
    else:
        if magnitude is None:
            raise RuntimeError("Effect magnitudes required when no chromatin-score cache exists")
        log_chrom = np.zeros(len(frame), dtype=np.float64)
        for feature in range(N_CHROMATIN):
            bg = magnitude[bg_ids, feature]
            ev = empirical_evalue(bg, magnitude[ids, feature])
            log_chrom += np.log(ev)
            if feature in (0, 299, 599, 918):
                print(f"{cohort}/{model}: empirical chromatin E-values {feature + 1}/{N_CHROMATIN}", flush=True)
        chrom_gm = np.exp(log_chrom / N_CHROMATIN)

    conservation_e = np.full((len(frame), len(CADD_COLUMNS)), np.nan, dtype=np.float64)


    polarity = np.ones(len(CADD_COLUMNS), dtype=np.float64)
    cadd = frame[list(CADD_COLUMNS)].to_numpy(dtype=np.float64) * polarity
    is_random = frame.label.eq("Random negative SNP").to_numpy()
    for idx, col in enumerate(CADD_COLUMNS):
        bg = cadd[is_random, idx]
        conservation_e[:, idx] = empirical_evalue(bg, cadd[:, idx])
    complete = np.isfinite(conservation_e).all(axis=1)
    conservation_gm = np.full(len(frame), np.nan, dtype=np.float64)
    conservation_gm[complete] = np.exp(np.log(conservation_e[complete]).mean(axis=1))
    combined = chrom_gm * conservation_gm

    out = frame[["chr", "pos", "ref", "alt", "label", "positive", "fold", "distance_bp", "_vid"]].copy()
    out = out.rename(columns={"_vid": "variant_id"})
    out["model"] = model
    for idx, col in enumerate(CADD_COLUMNS):
        out[col] = frame[col].to_numpy()
        out[f"E_{col}"] = conservation_e[:, idx]
    out["chromatin_geomean_E"] = chrom_gm
    out["conservation_geomean_E"] = conservation_gm
    out["functional_significance_score"] = combined
    out["all_four_conservation_scores_present"] = complete
    return out, bg_ids


def make_metrics(scored, cohort):
    rows = []
    is_positive = scored.positive.astype(bool).to_numpy()
    score_specs = [
        ("chromatin_geomean_E", "all cohort rows"),
        ("conservation_geomean_E", "four-conservation complete cases"),
        ("functional_significance_score", "four-conservation complete cases"),
    ]
    for group in NEGATIVE_GROUPS:
        neg = scored.label.eq(group).to_numpy()
        selected = is_positive | neg
        if not neg.any():
            continue
        for score, case_def in score_specs:
            values = scored[score].to_numpy(dtype=np.float64)
            valid = selected & np.isfinite(values)
            labels = is_positive[valid].astype(int)
            if len(np.unique(labels)) < 2:
                auc = np.nan
            else:
                auc = float(roc_auc_score(labels, -values[valid]))
            pos_mask, neg_mask = valid & is_positive, valid & neg
            rows.append({
                "cohort": cohort, "model": scored.model.iloc[0], "negative_group": group,
                "score": score, "case_definition": case_def,
                "n_positive": int(pos_mask.sum()), "n_negative": int(neg_mask.sum()),
                "positive_coverage": float(pos_mask.sum() / max(is_positive.sum(), 1)),
                "negative_coverage": float(neg_mask.sum() / max(neg.sum(), 1)),
                "auc": auc,
            })
    return rows


def main():
    sys.stdout.reconfigure(line_buffering=True)
    status = json.loads((HERE / "cadd_annotation_status.json").read_text())
    files, annotation, n_variants = load_annotations()
    if status["completed"] != status["total"] or len(files) != status["total"]:
        raise RuntimeError(f"Expected all CADD batches; status={status['completed']}/{status['total']}, files={len(files)}")
    batch_audit = json.loads((OUTPUTS / "cadd_downloaded_batch_audit.json").read_text())
    if batch_audit["downloaded_files"] != status["total"] or batch_audit["cross_batch_score_conflicts"]:
        raise RuntimeError("CADD batch audit is incomplete or has conflicting overlap scores")
    cadd_by_id = np.full((n_variants, len(CADD_COLUMNS)), np.nan, dtype=np.float32)
    ids = annotation.variant_id.to_numpy(dtype=np.int64)
    if len(np.unique(ids)) != len(ids):
        raise RuntimeError("CADD records map to duplicate experiment variant IDs")
    cadd_by_id[ids] = annotation[list(CADD_COLUMNS)].to_numpy(dtype=np.float32)
    universe_variants_without_cadd_record = n_variants - len(ids)
    variants_absent_from_cadd_vcf = n_variants - int(batch_audit["submitted_unique_variants"])
    submitted_variants_without_returned_cadd = int(batch_audit["submitted_unique_variants"] - len(annotation))
    returned_rows_missing_all_four = int(annotation[list(CADD_COLUMNS)].isna().all(axis=1).sum())
    cohorts = load_cohorts(cadd_by_id, n_variants)

    coverage_rows = []
    for cohort, frame in cohorts.items():
        full = frame[list(CADD_COLUMNS)].notna().all(axis=1)
        present = frame[list(CADD_COLUMNS)].notna().any(axis=1)
        for group, sub in frame.groupby("label", sort=False):
            ix = sub.index
            coverage_rows.append({
                "cohort": cohort, "label": group, "n_total": int(len(sub)),
                "n_any_downloaded_annotation": int(present.loc[ix].sum()),
                "n_complete_four_annotations": int(full.loc[ix].sum()),
                "any_coverage_fraction": float(present.loc[ix].mean()),
                "complete_coverage_fraction": float(full.loc[ix].mean()),
            })

    all_metrics = []
    provenance = {
        "analysis": "unsupervised functional significance using all returned CADD inclAnno batches",
        "as_of_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cadd_status_snapshot": {"completed": status["completed"], "total": status["total"],
                                 "last_checked_utc": status["last_checked_utc"]},
        "downloaded_batches": [p.name for p in files],
        "downloaded_unique_cadd_variants": int(len(annotation)),
        "analysis_effect_universe_variants": int(n_variants),
        "cadd_vcf_submitted_unique_variants": int(batch_audit["submitted_unique_variants"]),
        "analysis_universe_not_submitted_to_cadd": variants_absent_from_cadd_vcf,
        "submitted_variants_without_returned_cadd_record": submitted_variants_without_returned_cadd,
        "returned_variants_missing_all_four_cadd_fields": returned_rows_missing_all_four,
        "analysis_universe_variants_without_any_cadd_record": universe_variants_without_cadd_record,
        "analysis_universe_variants_without_any_usable_cadd_score": int(np.isnan(cadd_by_id).all(axis=1).sum()),
        "feature_source": str(FIGURE3.relative_to(HERE.parent.parent)),
        "feature_note": "Uses existing Figure3 local-inference effect matrices; no new model inference was run.",
        "models": {"pretrained": str(FIGURE3 / "pretrained_effects.npy"),
                   "ours": str(FIGURE3 / "ours_effects.npy")},
        "chromatin_effect_magnitude": "abs(P_ref-P_alt) * abs(logit(P_ref)-logit(P_alt)); 919 chromatin tracks",
        "background": "Unique Random negative SNP variants separately within each GRASP eQTL and GWAS cohort",
        "empirical_evalue": "(number of background scores strictly greater + 1)/(number of nonmissing background scores + 1); add-one smoothing for geometric means",
        "conservation_direction": "Scores use their reported direction; each E-value is the background proportion with a higher score, following DeepSEA.",
        "geometric_scores": "chromatin score = geometric mean of 919 E-values; conservation score = geometric mean of four E-values; combined score = their product",
        "metric_direction": "lower functional significance score indicates greater significance; AUC calculated using its negative",
        "limitations": [f"All {status['total']} CADD batches returned; {submitted_variants_without_returned_cadd:,} of {int(batch_audit['submitted_unique_variants']):,} submitted unique SNVs have no returned CADD record.",
                        f"The effect universe has {variants_absent_from_cadd_vcf:,} variants not included in the CADD VCF; these remain unannotated in this experiment.",
                        f"{returned_rows_missing_all_four:,} returned records lack all four required conservation values.",
                        "Individual conservation fields have missing values; conservation and combined AUCs use rows complete for all four fields.",
                        "These complete-case AUCs describe the returned cohort rows and are not unbiased estimates if missingness is systematic.",
                        "No classifier was trained; effect matrices were reused from Figure3 and no new inference was run."],
    }

    for model in ("pretrained", "ours"):
        cached = {}
        for cohort, frame in cohorts.items():
            prior_path = OUTPUTS / f"{cohort}_partial_unsupervised_{model}.tsv.gz"
            if not prior_path.exists():
                continue
            prior = pd.read_csv(prior_path, sep="\t", usecols=["variant_id", "chromatin_geomean_E"])
            if (len(prior) == len(frame) and
                    np.array_equal(prior.variant_id.to_numpy(dtype=np.int64), frame._vid.to_numpy(dtype=np.int64))):
                cached[cohort] = prior.chromatin_geomean_E.to_numpy(dtype=np.float64)
        magnitude = None
        if len(cached) != len(cohorts):
            effect_path = FIGURE3 / f"{model}_effects.npy"
            magnitude = prepare_effect_magnitude(effect_path, n_variants)
        for cohort, frame in cohorts.items():
            scored, bg_ids = score_cohort(frame, magnitude, model, cohort, cached.get(cohort))
            out_path = OUTPUTS / f"{cohort}_unsupervised_{model}.tsv.gz"
            scored.to_csv(out_path, sep="\t", index=False, compression="gzip", na_rep="NA")
            all_metrics.extend(make_metrics(scored, cohort))
            print(f"Wrote {out_path.name}; {len(scored):,} variants, {len(bg_ids):,} random-negative background variants", flush=True)
        if magnitude is not None:
            del magnitude

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(coverage_rows).to_csv(OUTPUTS / "annotation_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(all_metrics).to_csv(OUTPUTS / "unsupervised_auc.tsv", sep="\t", index=False, na_rep="NA")
    (OUTPUTS / "analysis_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Wrote annotation_coverage.tsv, unsupervised_auc.tsv, and analysis_manifest.json", flush=True)


if __name__ == "__main__":
    main()
