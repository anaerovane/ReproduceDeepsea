"""Draw separate Experiment 9 Supplementary Figure 6 panels for each model."""
from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
FIGURE3 = HERE.parent / "Figure3" / "outputs"
OUTPUTS = HERE / "outputs"
METRICS = OUTPUTS / "unsupervised_auc.tsv"
GROUPS = [
    "Random negative SNP",
    "31kbp negative SNP",
    "6.3kbp negative SNP",
    "710bp negative SNP",
    "360bp negative SNP",
]
COHORTS = [("eqtl", "GRASP eQTL (Noncoding)"), ("gwas", "GWAS Catalog (Noncoding)")]
STYLE = {"color": "#176b87", "linewidth": 2.2, "marker": None}


def draw_model(data: pd.DataFrame, model: str) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), sharey=False)
    for ax, (cohort, cohort_title) in zip(axes, COHORTS):
        sub = data[(data.cohort == cohort) &
                   (data.model == model) &
                   (data.score == "functional_significance_score")]
        by_group = sub.set_index("negative_group")
        aucs = [float(by_group.loc[group, "auc"]) for group in GROUPS]
        positive_n = int(by_group.loc[GROUPS[0], "n_positive"])
        summary = json.loads((FIGURE3 / f"{cohort}_pretrained_summary.json").read_text())
        distance_by_group = {row["group"]: row["mean_distance_bp"] for row in summary}
        tick_labels = ["All"] + [f"{distance_by_group[group]:,.0f}" for group in GROUPS[1:]]
        ax.plot(range(len(GROUPS)), aucs, **STYLE)
        ax.set_title(f"{cohort_title}\n(n = {positive_n:,} positive complete cases)", fontsize=10)
        ax.set_xticks(range(len(GROUPS)), tick_labels)
        ax.set_xlabel("Mean distance to nearest positive (bp)", fontsize=9)
        ax.set_ylabel("AUC", fontsize=9)
        ax.set_ylim(0.50, max(0.62, max(aucs) + 0.035))
        ax.set_xlim(-0.15, len(GROUPS) - 0.85)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.suptitle(f"DeepSEA functional significance score — {model}", fontsize=11, y=1.02)
    fig.tight_layout()
    out = OUTPUTS / f"suppfig6_{model}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    data = pd.read_csv(METRICS, sep="\t")
    for model in ("pretrained", "ours"):
        print(f"Wrote {draw_model(data, model)}")


if __name__ == "__main__":
    main()
