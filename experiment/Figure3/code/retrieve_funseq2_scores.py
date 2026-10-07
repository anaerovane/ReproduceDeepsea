"""Query the official hg19 FunSeq2 2.1.0 whole-genome track for cohort loci."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXPERIMENT_DIR = HERE.parent
OUTPUT_DIR = EXPERIMENT_DIR / 'outputs'
SOURCE = 'http://archive.gersteinlab.org/funseq2/hg19_wg_score.tsv.gz'
CHROM_ORDER = {f'chr{i}':i for i in range(1,23)} | {'chrX':23,'chrY':24}
RANDOM_NEGATIVE_FRACTION = 0.05
RANDOM_SEED = 42


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tabix',default='tabix')
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--track',default=SOURCE,
                        help='indexed FunSeq2 score track (local path or tabix-readable URL)')
    parser.add_argument('--output',type=Path,default=OUTPUT_DIR/'funseq2_v2.1.0_scores.tsv.gz')
    args=parser.parse_args()
    cohort=[]
    for task in ['eqtl','gwas']:
        frame=pd.read_csv(OUTPUT_DIR/f'{task}_pretrained_oof.tsv.gz',sep='\t')
        random_neg=frame[frame.label.eq('Random negative SNP')].sample(
            frac=RANDOM_NEGATIVE_FRACTION,random_state=RANDOM_SEED)
        other=frame[frame.positive|~frame.label.eq('Random negative SNP')]
        cohort.append(pd.concat([other,random_neg],ignore_index=True))
    variants=pd.concat(cohort,ignore_index=True)
    loci=variants[['chr','pos']].drop_duplicates()
    loci=loci[loci.chr.isin(CHROM_ORDER)].copy()
    loci['_chrom_order']=loci.chr.map(CHROM_ORDER)
    loci=loci.sort_values(['_chrom_order','pos'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    print(f'Querying {len(loci):,} unique cohort loci from FunSeq2 hg19 track; '
          f'All random negatives sampled at {RANDOM_NEGATIVE_FRACTION:.0%}',flush=True)
    def query_chrom(chrom, group, tempdir):
        region_path=Path(tempdir)/f'{chrom}.bed'
        with region_path.open('w') as f:
            for row in group.itertuples(index=False):
                f.write(f'{row.chr}\t{int(row.pos)-1}\t{int(row.pos)}\n')
        proc=subprocess.run([args.tabix,'-R',str(region_path),args.track],
                            capture_output=True,text=True)
        if proc.returncode:
            raise RuntimeError(f'{chrom}: {proc.stderr[-3000:]}')
        rows=[]
        for line in proc.stdout.splitlines():
            fields=line.split('\t')
            if len(fields)>=5 and not fields[0].startswith('#'):
                rows.append('\t'.join(fields[:5]))
        return chrom,rows

    total=0
    with tempfile.TemporaryDirectory(prefix='funseq2_regions_') as tempdir, gzip.open(args.output,'wt') as out:
        out.write('chr\tpos\tref\talt_set\tfunseq2_score\n')
        groups=list(loci.groupby('chr',sort=False))
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(query_chrom,chrom,group,tempdir) for chrom,group in groups]
            for future in as_completed(futures):
                chrom,rows=future.result()
                for row in rows:
                    out.write(row+'\n')
                total+=len(rows)
                out.flush()
                print(f'{chrom}: {len(rows):,} rows; total {total:,}',flush=True)
    count=total
    print(f'Wrote {count:,} FunSeq2 rows to {args.output}',flush=True)


if __name__=='__main__':
    main()
