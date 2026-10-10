import json
import time
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
MAPPING = DATA / 'grasp_variation_mappings.tsv.gz'
CLINVAR = DATA / 'clinvar_likely_pathogenic_GRCh37.tsv.gz'
OUTPUT = DATA / 'figure4d_positive_variants.tsv.gz'

def choose_grasp_alleles(mapping):
    parts = mapping.allele_string.fillna('').astype(str).str.split('/')
    keep = parts.map(len).ge(2)
    frame = mapping.loc[keep].copy()
    frame['ref'] = parts.loc[keep].map(lambda x: x[0].upper())
    frame['alt'] = parts.loc[keep].map(lambda x: [allele.upper() for allele in x[1:]])
    frame = frame.explode('alt', ignore_index=True)
    frame['minor_allele'] = frame.minor_allele.fillna('').astype(str).str.upper()
    frame['maf'] = pd.to_numeric(frame.maf, errors='coerce')
    return frame[frame.coordinate_match & frame.ref.isin(list('ACGT')) & frame.alt.isin(list('ACGT')) & frame.ref.ne(frame.alt) & frame.var_class.eq('SNP')].copy()

def main():
    mapping = pd.read_csv(MAPPING, sep='\t', dtype={'source_chrom': str, 'mapped_chrom': str})
    grasp = choose_grasp_alleles(mapping)
    grasp = grasp.rename(columns={'source_chrom': 'chrom', 'source_pos': 'pos', 'maf': 'dbsnp_maf'})
    grasp['positive'] = True
    grasp['source'] = 'GRASP 2.0 locus; dbSNP alternate-allele proxy'
    grasp['label'] = grasp.cohort.map({'eqtl': 'eQTL', 'meqtl': 'meQTL'})
    grasp['source_annotation'] = grasp['p_value'].astype(str)
    grasp = grasp[['cohort', 'label', 'chrom', 'pos', 'ref', 'alt', 'positive', 'source', 'rsid', 'source_annotation', 'dbsnp_maf']]
    grasp['pos'] = pd.to_numeric(grasp.pos, errors='coerce')
    clinvar = pd.read_csv(CLINVAR, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    clinvar = clinvar.rename(columns={'CLNSIG': 'source_annotation'})
    clinvar['cohort'] = 'clinvar'
    clinvar['label'] = 'ClinVar'
    clinvar['positive'] = True
    clinvar['source'] = 'ClinVar GRCh37 archive 2025-01-06'
    clinvar['rsid'] = clinvar.variant_id
    clinvar['dbsnp_maf'] = pd.NA
    clinvar['source_annotation'] = clinvar.source_annotation.astype(str)
    clinvar = clinvar[['cohort', 'label', 'chrom', 'pos', 'ref', 'alt', 'positive', 'source', 'rsid', 'source_annotation', 'dbsnp_maf']]
    output = pd.concat([grasp, clinvar], ignore_index=True)
    output['chrom'] = output.chrom.astype(str).str.replace('^chr', '', regex=True)
    output['pos'] = pd.to_numeric(output.pos, errors='coerce')
    output = output[output.chrom.isin({str(x) for x in range(1, 23)})].dropna(subset=['pos'])
    output['pos'] = output.pos.astype('int64')
    output = output.drop_duplicates(['cohort', 'chrom', 'pos', 'ref', 'alt'])
    output = output.sort_values(['cohort', 'chrom', 'pos', 'ref', 'alt'], kind='stable')
    output.insert(0, 'variant_id', range(len(output)))
    output.to_csv(OUTPUT, sep='\t', index=False, compression='gzip')
    summary = output.groupby(['cohort', 'label'], as_index=False).agg(positive_count=('variant_id', 'size'), chromosomes=('chrom', 'nunique'), dbsnp_maf_count=('dbsnp_maf', 'count'))
    summary.to_csv(DATA / 'positive_variant_audit.tsv', sep='\t', index=False)
    manifest = {'input_mapping': str(MAPPING), 'input_clinvar': str(CLINVAR), 'grasp_allele_rule': 'GRASP does not report the tested allele; retain each single-nucleotide alternate allele from coordinate-matched GRCh37 dbSNP mappings as a locus-level positive proxy', 'clinvar_filter': 'autosomal single-nucleotide alleles with CLNSIG containing Likely_pathogenic', 'positive_rows': len(output), 'output': str(OUTPUT), 'created_at_utc': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}
    (DATA / 'positive_variant_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.to_string(index=False), flush=True)
    print(f'Wrote {len(output):,} positives to {OUTPUT}', flush=True)
if __name__ == '__main__':
    main()
