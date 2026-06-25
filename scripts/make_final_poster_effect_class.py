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

REGIONS = {
    "5utr": {
        "label": "5\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
    },
    "tss": {
        "label": "Upstream Region",
        "summary": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
    },
    "3utr": {
        "label": "3\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
    },
}

EFFECT_MAP = {
    "qgrs_lost_all": "Complete G4 Loss",
    "lower_max_gscore": "G4 Weakened",
    "fewer_qgrs_lower_max_gscore": "G4 Weakened",
    "fewer_qgrs_same_max_gscore": "Motif Lost,\nScore Unchanged",
    "no_summary_change": "No Change",
    "more_qgrs_same_max_gscore": "Motif Gained,\nScore Unchanged",
    "higher_max_gscore": "G4 Strengthened",
}

CATEGORY_ORDER = [
    "Complete G4 Loss",
    "G4 Weakened",
    "Motif Lost,\nScore Unchanged",
    "No Change",
    "Motif Gained,\nScore Unchanged",
    "G4 Strengthened",
]

COLORS = {
    "Complete G4 Loss": "#b2182b",
    "G4 Weakened": "#ef8a62",
    "Motif Lost,\nScore Unchanged": "#fddbc7",
    "No Change": "#e0e0e0",
    "Motif Gained,\nScore Unchanged": "#d1e5f0",
    "G4 Strengthened": "#2166ac",
}


def load_effect_counts(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    df["effect_simple"] = df["EFFECT_CLASS"].map(EFFECT_MAP).fillna("Other")
    counts = df.groupby("effect_simple").size()
    return counts


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    all_counts = {}
    for key, meta in REGIONS.items():
        all_counts[meta["label"]] = load_effect_counts(meta["summary"])

    region_labels = list(all_counts.keys())
    x = np.arange(len(region_labels))
    width = 0.55

    fig, ax = plt.subplots(figsize=(10, 6))

    bottoms = np.zeros(len(region_labels))
    for cat in CATEGORY_ORDER:
        vals = []
        for rl in region_labels:
            total = all_counts[rl].sum()
            count = all_counts[rl].get(cat, 0)
            vals.append(100.0 * count / total)
        vals = np.array(vals)
        bars = ax.bar(x, vals, width, bottom=bottoms, label=cat,
                       color=COLORS[cat], edgecolor="white", linewidth=0.8)
        for i, v in enumerate(vals):
            if v >= 4:
                ax.text(x[i], bottoms[i] + v / 2, f"{v:.1f}%",
                        ha="center", va="center", fontsize=8.5, color="#222")
        bottoms += vals

    ax.set_xticks(x, region_labels, fontsize=11)
    ax.set_ylabel("% of G-related SNVs", fontsize=11)
    ax.set_ylim(0, 105)
    ax.set_title("Effect of Cancer-Associated SNVs on G-Quadruplex Motifs",
                  fontweight="bold", fontsize=13)
    ax.legend(loc="upper right", fontsize=9, frameon=True)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    outfile = OUTPUT_DIR / "combined_effect_class_breakdown.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
