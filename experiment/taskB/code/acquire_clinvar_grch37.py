import gzip
import hashlib
import json
import os
import time
from pathlib import Path
import pandas as pd
import pysam
import requests
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
URL = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2025/clinvar_20250106.vcf.gz'
VCF = DATA / 'clinvar_20250106_GRCh37.vcf.gz'
OUTPUT = DATA / 'clinvar_likely_pathogenic_GRCh37.tsv.gz'

def download():
    if VCF.is_file() and VCF.stat().st_size > 10000000:
        return
    VCF.parent.mkdir(parents=True, exist_ok=True)
    part = VCF.with_suffix(VCF.suffix + '.part')
    with requests.get(URL, stream=True, timeout=60) as response:
        response.raise_for_status()
        with part.open('wb') as handle:
            for chunk in response.iter_content(8 * 1024 * 1024):
                if chunk:
                    handle.write(chunk)
    if part.stat().st_size < 10000000:
        raise RuntimeError(f'ClinVar archive download is unexpectedly small: {part.stat().st_size}')
    part.replace(VCF)

def likely_pathogenic(value):
    values = value if isinstance(value, tuple) else (value,)
    text = '|'.join((str(item) for item in values)).casefold()
    text = text.replace(' ', '_')
    return 'likely_pathogenic' in text

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    download()
    records = []
    reader = pysam.VariantFile(str(VCF))
    for record in reader:
        chrom = str(record.contig).removeprefix('chr')
        if chrom not in {str(x) for x in range(1, 23)}:
            continue
        significance = record.info.get('CLNSIG')
        if not likely_pathogenic(significance):
            continue
        for alt in record.alts or ():
            if len(record.ref) != 1 or len(alt) != 1:
                continue
            if record.ref not in 'ACGT' or alt not in 'ACGT' or alt == record.ref:
                continue
            records.append({'chrom': chrom, 'pos': int(record.pos), 'ref': record.ref, 'alt': alt, 'variant_id': str(record.id or '.'), 'CLNSIG': '|'.join((str(x) for x in significance)) if isinstance(significance, tuple) else str(significance), 'CLNVC': str(record.info.get('CLNVC', '')), 'ALLELEID': str(record.info.get('ALLELEID', ''))})
    frame = pd.DataFrame.from_records(records)
    if frame.empty:
        raise RuntimeError('No likely pathogenic SNVs found in the GRCh37 ClinVar VCF')
    frame = frame.drop_duplicates(['chrom', 'pos', 'ref', 'alt'])
    frame.to_csv(OUTPUT, sep='\t', index=False, compression='gzip')
    digest = hashlib.sha256(VCF.read_bytes()).hexdigest()
    manifest = {'source_url': URL, 'snapshot': '2025-01-06', 'assembly': 'GRCh37', 'filter': 'autosomal single nucleotide variants with CLNSIG containing Likely_pathogenic', 'rows': len(frame), 'vcf_bytes': VCF.stat().st_size, 'vcf_sha256': digest, 'output': str(OUTPUT), 'updated_at_utc': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}
    (DATA / 'clinvar_grch37_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2), flush=True)
if __name__ == '__main__':
    main()
