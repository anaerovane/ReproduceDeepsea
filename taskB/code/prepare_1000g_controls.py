import concurrent.futures
import hashlib
import heapq
import json
import time
import threading
from bisect import bisect_left
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
POSITIVES = DATA / 'figure4d_positive_variants.tsv.gz'
SAMPLED_POSITIVES = DATA / 'figure4d_evaluation_positives.tsv.gz'
OUTPUT = DATA / 'figure4d_negative_variants.tsv.gz'
STATUS = DATA / 'negative_acquisition_status.json'
BASE_URL = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502'
SITES_VCF = f'{BASE_URL}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz'
COHORTS = ['eqtl', 'meqtl', 'clinvar']
FLANK_BP = 100000
MIN_MAF = 0.05
SEED = 20261010
SAMPLE_PER_COHORT = 10000
CHROM_WORKERS = 4
INTERVAL_WORKERS = 8
STOP_AFTER_CANDIDATES_PER_QUOTA = 20
_VCF_LOCAL = threading.local()

def vcf_url(chrom):
    return SITES_VCF

def merge_windows(positions, chrom_length=300000000):
    windows = sorted(((max(0, int(pos) - FLANK_BP - 1), min(chrom_length, int(pos) + FLANK_BP)) for pos in positions))
    merged = []
    for start, end in windows:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged

def nearest_distance(sorted_positions, pos):
    index = bisect_left(sorted_positions, pos)
    candidates = []
    if index < len(sorted_positions):
        candidates.append(abs(sorted_positions[index] - pos))
    if index:
        candidates.append(abs(sorted_positions[index - 1] - pos))
    return min(candidates) if candidates else None

def parse_af(info):
    for field in info.split(';'):
        if field.startswith('AF='):
            try:
                return [float(value) for value in field[3:].split(',')]
            except ValueError:
                return []
    return []

def variant_priority(cohort, chrom, pos, ref, alt):
    payload = f'{SEED}|{cohort}|{chrom}|{pos}|{ref}|{alt}'.encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), 'big')

def scan_interval(chrom, interval, positions, positive_position_sets, quotas):
    source = getattr(_VCF_LOCAL, 'source', None)
    if source is None:
        source = pysam.TabixFile(vcf_url(chrom))
        _VCF_LOCAL.source = source
    start, end = interval
    counts = {cohort: 0 for cohort in COHORTS}
    heaps = {cohort: [] for cohort in COHORTS}
    for line in source.fetch(chrom, start, end):
        fields = line.split('\t')
        pos = int(fields[1])
        ref = fields[3].upper()
        alts = fields[4].upper().split(',')
        if len(ref) != 1 or ref not in 'ACGT':
            continue
        afs = parse_af(fields[7])
        if len(afs) != len(alts):
            continue
        for alt, af in zip(alts, afs):
            if len(alt) != 1 or alt not in 'ACGT' or alt == ref:
                continue
            if not MIN_MAF < af < 1 - MIN_MAF:
                continue
            for cohort in COHORTS:
                quota = int(quotas.get(cohort, 0))
                if quota <= 0 or pos in positive_position_sets[cohort]:
                    continue
                distance = nearest_distance(positions[cohort], pos)
                if distance is None or distance > FLANK_BP:
                    continue
                counts[cohort] += 1
                row = {'cohort': cohort, 'label': {'eqtl': 'eQTL', 'meqtl': 'meQTL', 'clinvar': 'ClinVar'}[cohort], 'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'positive': False, 'source': '1000 Genomes Phase 3 v5c sites-only', 'rsid': fields[2], 'source_annotation': f'AF={af:.8g};MAF={min(af, 1 - af):.8g}', 'dbsnp_maf': min(af, 1 - af), 'distance_bp': int(distance)}
                priority = variant_priority(cohort, chrom, pos, ref, alt)
                heap = heaps[cohort]
                item = (-priority, pos, alt, row)
                if len(heap) < quota:
                    heapq.heappush(heap, item)
                elif item > heap[0]:
                    heapq.heapreplace(heap, item)
    return (counts, heaps)

def scan_chrom(chrom, cohort_positions, quotas, all_positive_positions):
    positions = {cohort: sorted(set((int(x) for x in cohort_positions[cohort].get(chrom, [])))) for cohort in COHORTS}
    all_positions = sorted(set().union(*(set(values) for values in positions.values())))
    if not all_positions:
        return {'chrom': chrom, 'negative_rows': [], 'candidates': {}, 'intervals': 0}
    reservoirs = {cohort: [] for cohort in COHORTS}
    seen = {cohort: 0 for cohort in COHORTS}
    positive_position_sets = {cohort: set(all_positive_positions[cohort].get(chrom, [])) for cohort in COHORTS}
    rng = {cohort: np.random.default_rng(SEED + int(chrom) * 17 + i) for i, cohort in enumerate(COHORTS)}
    intervals = merge_windows(all_positions)
    interval_rng = np.random.default_rng(SEED + int(chrom) * 101)
    interval_rng.shuffle(intervals)
    interval_pool = concurrent.futures.ThreadPoolExecutor(max_workers=INTERVAL_WORKERS)
    batch_size = INTERVAL_WORKERS * 4
    for offset in range(0, len(intervals), batch_size):
        batch = intervals[offset:offset + batch_size]
        futures = [interval_pool.submit(scan_interval, chrom, interval, positions, positive_position_sets, {c: int(quotas[c].get(chrom, 0)) for c in COHORTS}) for interval in batch]
        for future in concurrent.futures.as_completed(futures):
            counts, heaps = future.result()
            for cohort in COHORTS:
                seen[cohort] += counts[cohort]
                keep = int(quotas[cohort].get(chrom, 0))
                destination = reservoirs[cohort]
                for item in heaps[cohort]:
                    if len(destination) < keep:
                        heapq.heappush(destination, item)
                    elif item > destination[0]:
                        heapq.heapreplace(destination, item)
        print(json.dumps({'chrom': chrom, 'intervals_completed': min(offset + len(batch), len(intervals)), 'intervals_total': len(intervals), 'candidate_controls': seen}), flush=True)
        if all((int(quotas[c].get(chrom, 0)) == 0 or seen[c] >= int(quotas[c][chrom]) * STOP_AFTER_CANDIDATES_PER_QUOTA for c in COHORTS)):
            break
    interval_pool.shutdown(wait=False, cancel_futures=True)
    reservoirs = {cohort: [item[3] for item in heap] for cohort, heap in reservoirs.items()}
    candidates = {cohort: seen[cohort] for cohort in COHORTS}
    rows = [row for cohort in COHORTS for row in reservoirs[cohort]]
    return {'chrom': chrom, 'negative_rows': rows, 'candidates': candidates, 'intervals': len(intervals)}

def main():
    positives = pd.read_csv(POSITIVES, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    sampled = pd.concat([frame.sample(n=min(SAMPLE_PER_COHORT, len(frame)), random_state=SEED + i) for i, (_, frame) in enumerate(positives.groupby('cohort', sort=False))], ignore_index=True).sort_values(['cohort', 'chrom', 'pos', 'ref', 'alt'], kind='stable')
    sampled.to_csv(SAMPLED_POSITIVES, sep='\t', index=False, compression='gzip')
    cohort_positions = {}
    quotas = {}
    all_positive_positions = {}
    for cohort in COHORTS:
        frame = sampled[sampled.cohort.eq(cohort)]
        cohort_positions[cohort] = {chrom: group.pos.astype(int).tolist() for chrom, group in frame.groupby('chrom', sort=False)}
        quotas[cohort] = frame.chrom.value_counts().astype(int).to_dict()
        all_positive_positions[cohort] = {chrom: set(group.pos.astype(int)) for chrom, group in positives[positives.cohort.eq(cohort)].groupby('chrom', sort=False)}
    chroms = sorted(set().union(*(set(cohort_positions[c]) for c in COHORTS)), key=int)
    results = []
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=CHROM_WORKERS) as pool:
        future_map = {pool.submit(scan_chrom, chrom, cohort_positions, quotas, all_positive_positions): chrom for chrom in chroms}
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            completed += 1
            status = {'chromosomes_completed': completed, 'chromosomes_total': len(chroms), 'completed_chromosomes': sorted((x['chrom'] for x in results), key=int), 'intervals_processed': sum((x['intervals'] for x in results)), 'candidate_controls': {cohort: sum((x['candidates'][cohort] for x in results)) for cohort in COHORTS}, 'negative_rows_selected': {cohort: sum((sum((r['cohort'] == cohort for r in x['negative_rows'])) for x in results)) for cohort in COHORTS}, 'updated_at_utc': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}
            STATUS.write_text(json.dumps(status, indent=2) + '\n')
            print(json.dumps(status), flush=True)
    negative_rows = [row for result in results for row in result['negative_rows']]
    negatives = pd.DataFrame.from_records(negative_rows)
    if not negatives.empty:
        negatives = negatives.drop_duplicates(['cohort', 'chrom', 'pos', 'ref', 'alt'])
        negatives.insert(0, 'variant_id', range(len(negatives)))
    negatives.to_csv(OUTPUT, sep='\t', index=False, compression='gzip')
    rows = []
    for cohort in COHORTS:
        n_pos = len(sampled[sampled.cohort.eq(cohort)])
        cohort_neg = negatives[negatives.cohort.eq(cohort)] if not negatives.empty else negatives
        rows.append({'cohort': cohort, 'positive_count': n_pos, 'negative_count': len(cohort_neg), 'negative_candidates_seen': sum((x['candidates'][cohort] for x in results)), 'balanced': n_pos == len(cohort_neg)})
    audit = pd.DataFrame(rows)
    audit.to_csv(DATA / 'negative_variant_audit.tsv', sep='\t', index=False)
    manifest = {'source': '1000 Genomes Project Phase 3 integrated sites-only variant calls, v5c, GRCh37', 'source_url': SITES_VCF, 'minor_allele_frequency': f'>{MIN_MAF} in the INFO AF field', 'distance': f'within {FLANK_BP} bp of a nearest cohort-specific positive on the same chromosome', 'sampling': 'chromosome-stratified uniform hash-priority reservoir sample, one negative per positive per chromosome; query windows are shuffled and scanning stops after 20x the stratum quota to bound I/O', 'positive_sampling': f'random sample capped at {SAMPLE_PER_COHORT} positives per cohort before nearby control scanning; seed={SEED}', 'seed': SEED, 'positive_rows': len(sampled), 'positive_sample': str(SAMPLED_POSITIVES), 'rows': len(negatives), 'output': str(OUTPUT), 'audit': str(DATA / 'negative_variant_audit.tsv'), 'created_at_utc': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}
    (DATA / 'negative_variant_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(audit.to_string(index=False), flush=True)
    print(f'Wrote {len(negatives):,} negatives to {OUTPUT}', flush=True)
if __name__ == '__main__':
    main()
