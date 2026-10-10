import argparse
import concurrent.futures
import json
import time
from pathlib import Path
import pandas as pd
import requests
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = DATA
INPUTS = {'eqtl': DATA / 'eqtls_grasp2_p_lt_1e-12_hg19.tsv', 'meqtl': DATA / 'meqtls_grasp2_p_lt_1e-12_hg19.tsv'}
BASE = 'https://grch37.rest.ensembl.org/variation/human'

def fetch_batch(ids, timeout):
    response = requests.post(BASE, params={'content-type': 'application/json'}, json={'ids': ids}, timeout=timeout)
    if response.status_code == 429:
        time.sleep(float(response.headers.get('Retry-After', '2')))
        response = requests.post(BASE, params={'content-type': 'application/json'}, json={'ids': ids}, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    records = []
    for rsid in ids:
        item = payload.get(rsid, {})
        mappings = [mapping for mapping in item.get('mappings', []) if mapping.get('assembly_name') == 'GRCh37' and mapping.get('coord_system') == 'chromosome']
        for mapping in mappings:
            records.append({'rsid': rsid, 'mapped_chrom': str(mapping.get('seq_region_name', '')), 'mapped_pos': mapping.get('start'), 'allele_string': mapping.get('allele_string', ''), 'minor_allele': item.get('minor_allele', ''), 'maf': item.get('MAF'), 'var_class': item.get('var_class', ''), 'evidence': ','.join(item.get('evidence', [])), 'mapping_strand': mapping.get('strand')})
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch-size', type=int, default=200)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--timeout', type=int, default=120)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    variants = []
    for cohort, path in INPUTS.items():
        frame = pd.read_csv(path, sep='\t', dtype={'chrom_hg19': str, 'rsid_dbsnp134': str})
        frame['rsid'] = 'rs' + frame.rsid_dbsnp134.astype(str).str.replace('^rs', '', regex=True)
        frame['source_chrom'] = frame.chrom_hg19.astype(str).str.replace('^chr', '', regex=True)
        frame['source_pos'] = pd.to_numeric(frame.pos_hg19, errors='coerce')
        frame['cohort'] = cohort
        variants.append(frame)
    source = pd.concat(variants, ignore_index=True)
    batches = [source.rsid.iloc[start:start + args.batch_size].drop_duplicates().tolist() for start in range(0, len(source), args.batch_size)]
    unique_batches = []
    seen = set()
    for batch in batches:
        new_batch = [rsid for rsid in batch if rsid not in seen]
        seen.update(new_batch)
        if new_batch:
            unique_batches.append(new_batch)
    all_records = []
    failed = []
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch_batch, batch, args.timeout): batch for batch in unique_batches}
        for future in concurrent.futures.as_completed(futures):
            batch = futures[future]
            try:
                all_records.extend(future.result())
            except Exception as exc:
                failed.append({'ids': batch, 'error': f'{type(exc).__name__}: {exc}'})
            completed += 1
            if completed % 10 == 0 or completed == len(futures):
                status = {'unique_rsids': len(seen), 'batches_completed': completed, 'batches_total': len(futures), 'mapping_rows_received': len(all_records), 'failed_batches': len(failed), 'updated_at': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}
                (OUT / 'grasp_mapping_status.json').write_text(json.dumps(status, indent=2) + '\n')
                print(json.dumps(status), flush=True)
    mapping = pd.DataFrame.from_records(all_records)
    if mapping.empty:
        raise RuntimeError('Ensembl returned no GRCh37 variation mappings')
    source = source.merge(mapping, on='rsid', how='left', validate='many_to_many')
    source['coordinate_match'] = source.source_chrom.eq(source.mapped_chrom) & pd.to_numeric(source.source_pos, errors='coerce').eq(pd.to_numeric(source.mapped_pos, errors='coerce'))
    source.to_csv(OUT / 'grasp_variation_mappings.tsv.gz', sep='\t', index=False, compression='gzip')
    mapping.to_csv(OUT / 'grch37_variation_lookup.tsv.gz', sep='\t', index=False, compression='gzip')
    (OUT / 'grasp_mapping_failures.json').write_text(json.dumps(failed, indent=2) + '\n')
    summary = source.groupby('cohort', dropna=False).agg(source_rows=('rsid', 'size'), mapped_rows=('mapped_pos', 'count'), coordinate_matches=('coordinate_match', 'sum'), minor_allele_known=('minor_allele', lambda x: x.notna().sum()), maf_known=('maf', lambda x: x.notna().sum())).reset_index()
    summary.to_csv(OUT / 'grasp_mapping_audit.tsv', sep='\t', index=False)
    print(summary.to_string(index=False), flush=True)
    print(f'Wrote cohort mapping audit to {OUT}', flush=True)
if __name__ == '__main__':
    main()
