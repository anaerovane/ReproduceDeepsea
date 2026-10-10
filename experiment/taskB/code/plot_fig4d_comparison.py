from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc, roc_curve
ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / 'outputs'
DATA = ROOT / 'data'
FUNSEQ = OUTPUTS / 'fig4d_funseq2_scores.tsv.gz'
PLOT_DIR = ROOT / 'figures'
KEYS = ['cohort', 'chrom', 'pos', 'ref', 'alt']
METHODS = [('NT-v2 50M fine-tuned', 'nt_finetuned_score', '#e76f51'), ('DeepSEA pretrained', 'deepsea_pretrained_sum_abs_delta', '#277da1'), ('DeepSEA ours', 'deepsea_ours_sum_abs_delta', '#43aa8b'), ('FunSeq2 v2.1.0', 'funseq2_score', '#8e6cbb')]
COHORTS = [('eqtl', 'eQTL'), ('meqtl', 'meQTL'), ('clinvar', 'ClinVar')]

def read_funseq_scores():
    scores = pd.read_csv(FUNSEQ, sep='\t', dtype={'chrom': str, 'chr': str, 'ref': str, 'alt_set': str})
    if 'chrom' not in scores and 'chr' in scores:
        scores = scores.rename(columns={'chr': 'chrom'})
    scores['chrom'] = scores['chrom'].str.removeprefix('chr')
    scores['pos'] = pd.to_numeric(scores['pos'], errors='coerce')
    scores['funseq2_score'] = pd.to_numeric(scores['funseq2_score'], errors='coerce')
    scores = scores.dropna(subset=['pos', 'funseq2_score']).copy()
    scores['pos'] = scores['pos'].astype('int64')
    exact = scores[scores.alt_set.ne('.')][['chrom', 'pos', 'ref', 'alt_set', 'funseq2_score']]
    exact = exact.rename(columns={'alt_set': 'alt'})
    generic = scores[scores.alt_set.eq('.')][['chrom', 'pos', 'ref', 'funseq2_score']]
    generic = generic.drop_duplicates(['chrom', 'pos', 'ref'])
    return (exact, generic)

def read_joined_scores():
    nt = pd.read_csv(OUTPUTS / 'nt50m_fig4d_oof.tsv.gz', sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    pretrained = pd.read_csv(OUTPUTS / 'deepsea_pretrained_scores.tsv.gz', sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    ours = pd.read_csv(OUTPUTS / 'deepsea_ours_scores.tsv.gz', sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    joined = nt.merge(pretrained[KEYS + ['deepsea_pretrained_sum_abs_delta']], on=KEYS, how='left', validate='one_to_one')
    joined = joined.merge(ours[KEYS + ['deepsea_ours_sum_abs_delta']], on=KEYS, how='left', validate='one_to_one')
    exact, generic = read_funseq_scores()
    joined = joined.merge(exact, on=['chrom', 'pos', 'ref', 'alt'], how='left', validate='many_to_one')
    fallback = joined[['chrom', 'pos', 'ref']].merge(generic, on=['chrom', 'pos', 'ref'], how='left', validate='many_to_one')['funseq2_score']
    joined['funseq2_score'] = joined['funseq2_score'].fillna(fallback)
    return joined

def make_figure(joined):
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    colors = {label: color for label, _, color in METHODS}
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.3), sharex=True, sharey=True)
    rows = []
    for ax, (cohort, title) in zip(axes, COHORTS):
        frame = joined[joined.cohort.eq(cohort)]
        common = frame.dropna(subset=[column for _, column, _ in METHODS]).copy()
        unavailable = []
        for label, column, color in METHODS:
            valid = common.dropna(subset=[column])
            positives = int(valid.positive.astype(bool).sum())
            negatives = int((~valid.positive.astype(bool)).sum())
            if valid.positive.astype(bool).nunique() < 2:
                unavailable.append(f'{label} (P{positives:,}/N{negatives:,})')
                rows.append({'cohort': cohort, 'method': label, 'positive_count': positives, 'negative_count': negatives, 'scored_count': len(valid), 'cohort_candidate_count': len(frame), 'all_method_common_count': len(common), 'auc': None, 'status': 'AUC unavailable: both classes are required'})
                continue
            fpr, tpr, _ = roc_curve(valid.positive.astype(int), valid[column].astype(float))
            value = float(auc(fpr, tpr))
            ax.plot(fpr, tpr, color=color, linewidth=2.0, label=f'{label}  AUC {value:.3f}  (P{positives:,}/N{negatives:,})')
            rows.append({'cohort': cohort, 'method': label, 'positive_count': positives, 'negative_count': negatives, 'scored_count': len(valid), 'cohort_candidate_count': len(frame), 'all_method_common_count': len(common), 'auc': value, 'status': 'ROC on identical variants scored by all four methods'})
        if unavailable:
            ax.text(0.04, 0.96, 'AUC unavailable: ' + ', '.join(unavailable), transform=ax.transAxes, color='#333333', fontsize=7.5, va='top')
        ax.plot([0, 1], [0, 1], color='#777777', linestyle='--', linewidth=1)
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.01)
        ax.set_xlabel('False positive rate')
        ax.grid(alpha=0.18, linewidth=0.6)
        ax.legend(loc='lower right', fontsize=7.2, frameon=True)
    axes[0].set_ylabel('Sensitivity')
    fig.suptitle('Figure 4d method comparison', y=1.04, fontsize=14)
    fig.text(0.5, 0.99, 'Same scored variants per panel; FunSeq2 uses its hg19 noncoding score track', ha='center', va='top', fontsize=9, color='#444444')
    fig.tight_layout()
    fig.savefig(PLOT_DIR / 'figure4d_method_comparison.png', dpi=220, bbox_inches='tight')
    fig.savefig(PLOT_DIR / 'figure4d_method_comparison.pdf', bbox_inches='tight')
    pd.DataFrame(rows).to_csv(OUTPUTS / 'figure4d_method_metrics.tsv', sep='\t', index=False)
    plt.close(fig)
if __name__ == '__main__':
    make_figure(read_joined_scores())
