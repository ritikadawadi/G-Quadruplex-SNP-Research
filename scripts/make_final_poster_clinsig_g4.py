#!/usr/bin/env python3
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

SUMMARIES = [
    PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
    PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
    PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
]

SIG_MAP = {
    "Pathogenic": "Pathogenic /\nLikely Pathogenic",
    "Likely_pathogenic": "Pathogenic /\nLikely Pathogenic",
    "Uncertain_significance": "Uncertain\nSignificance",
    "Conflicting_classifications_of_pathogenicity": "Conflicting",
    "Likely_benign": "Benign /\nLikely Benign",
    "Benign": "Benign /\nLikely Benign",
    "Benign/Likely_benign": "Benign /\nLikely Benign",
}

CAT_ORDER = [
    "Pathogenic /\nLikely Pathogenic",
    "Uncertain\nSignificance",
    "Conflicting",
    "Benign /\nLikely Benign",
]

COLORS = {
    "Pathogenic /\nLikely Pathogenic": "#d7301f",
    "Uncertain\nSignificance": "#fc8d59",
    "Conflicting": "#fee090",
    "Benign /\nLikely Benign": "#4575b4",
}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    frames = []
    for p in SUMMARIES:
        frames.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(frames, ignore_index=True)

    df["sig_group"] = df["CLNSIG_PRIMARY"].map(SIG_MAP)
    df = df.dropna(subset=["sig_group"])
    df["changed"] = df["DELTA_MAX_GSCORE"] != 0

    stats = df.groupby("sig_group").agg(
        n=("changed", "count"),
        n_changed=("changed", "sum"),
        mean_abs_delta=("DELTA_MAX_GSCORE", lambda s: s[s != 0].abs().mean() if (s != 0).any() else 0),
    ).reindex(CAT_ORDER).dropna()

    stats["pct_changed"] = 100.0 * stats["n_changed"] / stats["n"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    x = np.arange(len(stats))
    colors = [COLORS[c] for c in stats.index]
    width = 0.55

    ax1.bar(x, stats["pct_changed"], width, color=colors, edgecolor="white", linewidth=0.8)
    for i, (cat, row) in enumerate(stats.iterrows()):
        ax1.text(i, row["pct_changed"] + 1.5, f"n={int(row['n'])}",
                 ha="center", fontsize=8.5, color="#555")
    ax1.set_xticks(x, stats.index, fontsize=9)
    ax1.set_ylabel("% of variants with\nmax G-score change", fontsize=10)
    ax1.set_title("How often does G-score change?", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, max(stats["pct_changed"]) * 1.3)
    ax1.grid(axis="y", alpha=0.25)

    ax2.bar(x, stats["mean_abs_delta"], width, color=colors, edgecolor="white", linewidth=0.8)
    for i, (cat, row) in enumerate(stats.iterrows()):
        n_ch = int(row["n_changed"])
        ax2.text(i, row["mean_abs_delta"] + 0.3, f"n={n_ch}",
                 ha="center", fontsize=8.5, color="#555")
    ax2.set_xticks(x, stats.index, fontsize=9)
    ax2.set_ylabel("Mean |Δ max G-score|\n(among changed variants)", fontsize=10)
    ax2.set_title("How severe is the disruption?", fontsize=11, fontweight="bold")
    ax2.grid(axis="y", alpha=0.25)

    fig.suptitle("Clinical Significance and G-Quadruplex Disruption",
                 fontsize=13, fontweight="bold", y=1.01)

    fig.text(0.5, -0.02,
             "All three regions (5\u2032UTR, upstream, 3\u2032UTR) combined. "
             "Clinical significance from ClinVar CLNSIG.",
             ha="center", fontsize=9, color="#666")

    fig.tight_layout()
    outfile = OUTPUT_DIR / "combined_clinsig_g4_disruption.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
