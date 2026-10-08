"""Plot AUCs for each Supplementary Figure 6 score across negative cohorts."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"
METRICS = OUTPUTS / "unsupervised_auc.tsv"
OUTPUT = OUTPUTS / "unsupervised_auc.png"
GROUPS = ["Random negative SNP", "31kbp negative SNP", "6.3kbp negative SNP", "710bp negative SNP", "360bp negative SNP"]
SCORES = [
    ("chromatin_geomean_E", "Chromatin"),
    ("conservation_geomean_E", "Conservation"),
    ("functional_significance_score", "Combined"),
]


def main():
    data = pd.read_csv(METRICS, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    x = range(len(GROUPS))
    styles = {"pretrained": ("#2463a6", "o", "-"), "ours": ("#d16b27", "s", "--")}
    for ax, cohort, title in zip(axes, ("eqtl", "gwas"), ("GRASP eQTL", "GWAS Catalog")):
        for column, label in SCORES:
            subset = data[(data.cohort == cohort) & data.negative_group.isin(GROUPS) & (data.score == column)]
            for model, (color, marker, linestyle) in styles.items():
                values = subset[subset.model == model].set_index("negative_group").reindex(GROUPS).auc
                ax.plot(x, values, color=color, marker=marker, linestyle=linestyle, linewidth=1.8,
                        markersize=5, label=f"{label} · {model}")
        ax.axhline(0.5, color="#777777", linewidth=1, linestyle=":")
        ax.set_title(title)
        ax.set_xticks(list(x), ["Random", "31kbp", "6.3kbp", "710bp", "360bp"], rotation=25, ha="right")
        ax.set_ylim(0.48, 0.66)
        ax.grid(axis="y", alpha=0.22)
        ax.set_xlabel("Negative variant group")
    axes[0].set_ylabel("AUC (positive vs. negative variants)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Supplementary Figure 6: unsupervised functional-significance AUC", y=1.02)
    fig.tight_layout(rect=(0, 0.12, 1, 0.98))
    fig.savefig(OUTPUT, dpi=180, bbox_inches="tight")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
