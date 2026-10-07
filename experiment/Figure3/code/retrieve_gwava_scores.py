#!/usr/bin/env python3
"""Retrieve indexed GRCh37 GWAVA scores for Figure3 cohort loci."""
import argparse
import gzip
import subprocess
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXPERIMENT_DIR = HERE.parent
OUT = EXPERIMENT_DIR / 'outputs'
RANDOM_FRACTION = 0.05
SEED = 42
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tabix', default='tabix')
    p.add_argument('--track', type=Path, default=OUT/'GRCh37_gwava.bed.gz')
    p.add_argument('--output', type=Path, default=OUT/'gwava_GRCh37_scores.tsv.gz')
    args = p.parse_args()

    target = []
    for task in ['eqtl', 'gwas']:
        frame = pd.read_csv(OUT/f'{task}_pretrained_oof.tsv.gz', sep='\t')
        random_neg = frame[frame.label.eq('Random negative SNP')].sample(
            frac=RANDOM_FRACTION, random_state=SEED)
        retained = frame[frame.positive | ~frame.label.eq('Random negative SNP')]
        target.append(pd.concat([retained, random_neg], ignore_index=True))
    loci = pd.concat(target, ignore_index=True)[['chr', 'pos']].drop_duplicates()
    loci = loci[loci.chr.isin(CHROMS)].copy()
    loci.pos = loci.pos.astype(int)
    print(f'Querying {len(loci):,} unique GRCh37 loci', flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    total_rows = 0
    with gzip.open(args.output, 'wt') as out:


        out.write('chr\tpos\tGWAVA_Region\tGWAVA_TSS\tGWAVA_Unmatched\ttrack_rows\tambiguous\n')
        for chrom in CHROMS:
            points = loci.loc[loci.chr.eq(chrom), 'pos'].drop_duplicates().sort_values()
            if points.empty:
                continue
            region_file = args.output.parent/f'.gwava_{chrom}.bed'
            with region_file.open('w') as bed:
                for pos in points:
                    bed.write(f'{chrom}\t{pos-1}\t{pos}\n')
            proc = subprocess.run([args.tabix, '-R', str(region_file), str(args.track)],
                                  capture_output=True, text=True)
            region_file.unlink(missing_ok=True)
            if proc.returncode:
                raise RuntimeError(f'tabix failed on {chrom}: {proc.stderr[-2000:]}')
            by_pos = {}
            for line in proc.stdout.splitlines():
                fields = line.split('\t')
                if len(fields) < 7:
                    continue
                c, start = fields[0], int(fields[1])
                key = (c, start+1)
                vals = tuple(float(x) for x in fields[4:7])
                by_pos.setdefault(key, set()).add(vals)
            for (c, pos), vectors in by_pos.items():


                ambiguous = len(vectors) > 1
                vals = next(iter(vectors)) if not ambiguous else (float('nan'),)*3
                out.write(f'{c}\t{pos}\t{vals[0]}\t{vals[1]}\t{vals[2]}\t{len(vectors)}\t{int(ambiguous)}\n')
                total_rows += 1
            out.flush()
            print(f'{chrom}: scored {len(by_pos):,}/{len(points):,} loci', flush=True)
    print(f'Wrote {total_rows:,} coordinate scores to {args.output}', flush=True)


if __name__ == '__main__':
    main()
