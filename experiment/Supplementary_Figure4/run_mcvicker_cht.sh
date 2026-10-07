set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WASP="${WASP_CHT_DIR:-}"
OUT="$HERE/results"
mkdir -p "$OUT"

if [ -z "$WASP" ] || [ ! -f "$WASP/combined_test.py" ]; then
  echo "Set WASP_CHT_DIR to the directory containing WASP combined_test.py" >&2
  exit 2
fi

for MARK in H3K4me3 H3K27ac; do
  LIST="$OUT/${MARK}_input_files.txt"
  find "$HERE/inputs" -maxdepth 1 -type f -name "*_$(printf '%s' "$MARK")_*_read_counts.txt.gz" -print | sort > "$LIST"
  if [ "$(wc -l < "$LIST")" -ne 10 ]; then
    echo "Expected 10 GEO donor files for $MARK; found $(wc -l < "$LIST")" >&2
    exit 2
  fi
done

if [ ! -s "$OUT/H3K4me3_as_dispersion.txt" ] || [ ! -s "$OUT/H3K27ac_as_dispersion.txt" ]; then
  python "$HERE/fit_as_dispersions_parallel.py"
fi
for MARK in H3K4me3 H3K27ac; do
  LIST="$OUT/${MARK}_input_files.txt"
  python "$WASP/fit_bnb_coefficients.py" --min_counts 50 --min_as_counts 10 --seed 3547 "$LIST" "$OUT/${MARK}_bnb_dispersion.txt"
  python "$WASP/combined_test.py" --min_as_counts 10 \
    --bnb_disp "$OUT/${MARK}_bnb_dispersion.txt" \
    --as_disp "$OUT/${MARK}_as_dispersion.txt" \
    "$LIST" "$OUT/${MARK}_cht_results.tsv"
done
python "$HERE/summarize_cht_qtls.py"
python "$HERE/reproduce_histone_qtl.py"
