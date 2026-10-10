import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd
import pysam
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUTPUTS = ROOT / 'outputs'
POSITIVES = DATA / 'figure4d_evaluation_positives.tsv.gz'
NEGATIVES = DATA / 'figure4d_negative_variants.tsv.gz'
TRACK = Path(os.environ.get('FUNSEQ2_TRACK', '/tmp/funseq2_fig4/hg19_wg_score.tsv.gz'))
OUTPUT = OUTPUTS / 'fig4d_funseq2_scores.tsv.gz'
AUDIT = OUTPUTS / 'fig4d_funseq2_coverage.tsv'
CLINVAR_VCF = DATA / 'clinvar_20250106_GRCh37.vcf.gz'
CONSEQUENCES = OUTPUTS / 'fig4d_clinvar_consequences.tsv.gz'
MANIFEST = OUTPUTS / 'fig4d_funseq2_manifest.json'
SOURCE_URL = 'http://archive.gersteinlab.org/funseq2/hg19_wg_score.tsv.gz'

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def write_clinvar_consequences(variants):
    positives = variants[(variants.cohort == 'clinvar') & variants.positive]
    wanted = {(row.chrom, int(row.pos), row.ref, row.alt): row.variant_id for row in positives.itertuples()}
    consequence_by_id = {}
    pysam.set_verbosity(0)
    vcf = pysam.VariantFile(str(CLINVAR_VCF))
    for record in vcf:
        chrom = record.contig.removeprefix('chr')
        for alt in record.alts or ():
            variant_id = wanted.get((chrom, int(record.pos), record.ref, alt))
            if variant_id is not None:
                consequence_by_id[variant_id] = ';'.join(record.info.get('MC') or ())
    annotated = positives[['variant_id', 'chrom', 'pos', 'ref', 'alt']].copy()
    annotated['cohort'] = 'clinvar'
    annotated['molecular_consequence'] = annotated.variant_id.map(consequence_by_id).fillna('not_in_vcf')
    coding_terms = ('missense_variant', 'nonsense', 'splice_donor_variant', 'splice_acceptor_variant', 'synonymous_variant', 'stop_lost', 'initiator_codon_variant', 'coding_sequence_variant', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion', 'transcript_ablation')
    text = annotated.molecular_consequence.str.lower()
    annotated['region_class'] = 'noncoding_or_other'
    annotated.loc[text.str.contains('not_in_vcf'), 'region_class'] = 'not_in_vcf'
    annotated.loc[text.apply(lambda value: any((term in value for term in coding_terms))), 'region_class'] = 'coding'
    annotated.to_csv(CONSEQUENCES, sep='\t', index=False, compression='gzip')
    return annotated

def main():
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    variants = pd.concat([pd.read_csv(path, sep='\t', dtype={'chrom': str}) for path in [POSITIVES, NEGATIVES]], ignore_index=True)
    variants['chrom'] = variants['chrom'].str.removeprefix('chr')
    consequence_table = write_clinvar_consequences(variants)
    positions = variants[['chrom', 'pos']].drop_duplicates().copy()
    positions['pos'] = positions['pos'].astype(int)
    rows = []
    track = pysam.TabixFile(str(TRACK))
    for chrom, group in positions.groupby('chrom', sort=False):
        contig = f'chr{chrom}'
        if contig not in track.contigs:
            continue
        for pos in group.pos:
            try:
                for line in track.fetch(contig, int(pos) - 1, int(pos)):
                    fields = line.split('\t')
                    if len(fields) < 5 or int(fields[1]) != int(pos):
                        continue
                    rows.append({'chrom': chrom, 'pos': int(fields[1]), 'ref': fields[2], 'alt_set': fields[3], 'funseq2_score': float(fields[4])})
            except (ValueError, OSError):
                continue
    scores = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt_set', 'funseq2_score'])
    DATA.mkdir(parents=True, exist_ok=True)
    scores.to_csv(OUTPUT, sep='\t', index=False, compression='gzip')
    exact = scores[scores.alt_set.ne('.')].rename(columns={'alt_set': 'alt'})
    generic = scores[scores.alt_set.eq('.')].drop_duplicates(['chrom', 'pos', 'ref'])
    joined = variants.merge(exact[['chrom', 'pos', 'ref', 'alt', 'funseq2_score']], on=['chrom', 'pos', 'ref', 'alt'], how='left', validate='many_to_one')
    fallback = variants.merge(generic[['chrom', 'pos', 'ref', 'funseq2_score']], on=['chrom', 'pos', 'ref'], how='left', validate='many_to_one')['funseq2_score']
    joined['score'] = joined['funseq2_score'].fillna(fallback)
    joined['scored'] = joined.score.notna()
    joined = joined.merge(consequence_table[['cohort', 'variant_id', 'region_class']], on=['cohort', 'variant_id'], how='left', validate='many_to_one')
    joined['region_class'] = joined.region_class.fillna('negative_control')
    audit = joined.groupby(['cohort', 'positive', 'region_class'], as_index=False).agg(total=('scored', 'size'), scored=('scored', 'sum'))
    audit['unscored'] = audit.total - audit.scored
    audit['coverage'] = audit.scored / audit.total
    audit.to_csv(AUDIT, sep='\t', index=False)
    MANIFEST.write_text(json.dumps({'method': 'FunSeq2 v2.1.0 hg19 precomputed noncoding score track', 'source_url': SOURCE_URL, 'source_track_bytes': TRACK.stat().st_size, 'source_index_url': SOURCE_URL + '.tbi', 'source_score_granularity': "the published whole-genome track may use generic alternate '.' rows; these are matched by chromosome, position, and reference allele", 'coordinate_matching': "exact hg19 chromosome, 1-based position, and reference; exact alternate when present, generic-alt fallback for track rows marked '.'", 'comparability': 'ROCs use variant rows scored by all four comparison methods', 'inputs': {str(POSITIVES.name): sha256(POSITIVES), str(NEGATIVES.name): sha256(NEGATIVES), str(CLINVAR_VCF.name): sha256(CLINVAR_VCF)}, 'outputs': {str(OUTPUT.name): sha256(OUTPUT), str(AUDIT.name): sha256(AUDIT), str(CONSEQUENCES.name): sha256(CONSEQUENCES)}, 'extraction_code_sha256': sha256(Path(__file__).resolve()), 'created_at_utc': datetime.now(timezone.utc).isoformat()}, indent=2) + '\n')
    print(audit.to_string(index=False))
    print('ClinVar positive consequence classes:')
    print(consequence_table.region_class.value_counts().to_string())
    print(f'Wrote {OUTPUT} ({len(scores):,} rows)')
if __name__ == '__main__':
    main()
