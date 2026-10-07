#!/usr/bin/env python3
"""Audit downloaded CADD annotations against submitted variants and each other."""
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
OUTPUTS = HERE / "outputs"
VCF = INPUTS / "eqtl_gwas_variants_for_cadd_inclAnno.vcf.gz"
OUT = OUTPUTS / "cadd_downloaded_batch_audit.json"
REQUIRED = ("priPhCons", "priPhyloP", "GerpN", "GerpS")


def key(chrom, pos, ref, alt):
    return (chrom.removeprefix("chr"), pos, ref.upper(), alt.upper())


def read_vcf():
    variants = set()
    with gzip.open(VCF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _vid, ref, alt, *_ = line.rstrip("\n").split("\t")
            variants.add(key(chrom, pos, ref, alt))
    return variants


def read_batch(path):
    rows = {}
    duplicate_rows = 0
    missing_values = {col: 0 for col in REQUIRED}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("##")), delimiter="\t")
        missing_cols = [col for col in REQUIRED if col not in (reader.fieldnames or [])]
        if missing_cols:
            raise ValueError(f"{path.name}: missing columns {missing_cols}")
        for row in reader:
            k = key(row["#Chrom"], row["Pos"], row["Ref"], row["Alt"])
            scores = tuple(row[col] for col in REQUIRED)
            for col, value in zip(REQUIRED, scores):
                if value in ("", "NA", "."):
                    missing_values[col] += 1
            if k in rows:
                duplicate_rows += 1
                if rows[k] != scores:
                    raise ValueError(f"Conflicting duplicate variant in {path.name}: {k}")
            rows[k] = scores
    return rows, duplicate_rows, missing_values


def main():
    submitted = read_vcf()
    files = sorted(OUTPUTS.glob("cadd_inclAnno_batch_*.tsv.gz"))
    batches = {}
    result_files = {}
    for path in files:
        rows, within_duplicates, missing = read_batch(path)
        keys = set(rows)
        unmatched = keys - submitted
        batches[path.name] = rows
        result_files[path.name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rows": sum(1 for _ in rows) + within_duplicates,
            "unique_variants": len(keys),
            "within_batch_duplicate_rows": within_duplicates,
            "required_columns_present": True,
            "rows_matching_submitted_vcf": len(keys) - len(unmatched),
            "unmatched_variants": len(unmatched),
            "missing_values": missing,
        }

    pairwise = {}
    for left, right in itertools.combinations(files, 2):
        a, b = batches[left.name], batches[right.name]
        shared = set(a) & set(b)
        different = sum(a[k] != b[k] for k in shared)
        pairwise[f"{left.name} vs {right.name}"] = {
            "shared_variants": len(shared),
            "different_required_score_tuples": different,
            "unique_to_left": len(set(a) - set(b)),
            "unique_to_right": len(set(b) - set(a)),
        }

    merged = {}
    conflicts = []
    for name, rows in batches.items():
        for k, scores in rows.items():
            if k in merged and merged[k] != scores:
                conflicts.append(k)
            else:
                merged[k] = scores
    audit = {
        "source": "CADD GRCh37-v1.4 inclAnno output; audited against submitted VCF by CHROM/POS/REF/ALT",
        "submitted_unique_variants": len(submitted),
        "downloaded_files": len(files),
        "downloaded_batch_files": result_files,
        "pairwise_overlaps": pairwise,
        "raw_rows_across_downloads": sum(v["rows"] for v in result_files.values()),
        "unique_variants_across_downloads": len(merged),
        "duplicate_rows_across_batches": sum(v["rows"] for v in result_files.values()) - len(merged),
        "cross_batch_score_conflicts": len(conflicts),
        "conclusion": "All downloaded annotations are keyed and deduplicated by chromosome, position, REF, and ALT; see per-file and pairwise audits.",
    }
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
