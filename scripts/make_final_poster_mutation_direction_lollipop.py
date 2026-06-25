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

CLR_FROM = "#d7301f"
CLR_TO = "#2171b5"


def load_direction_stats(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t")
    out = {}
    for cls in ["from_G", "to_G"]:
        sub = df[df["G_CHANGE_CLASS"] == cls]
        n = len(sub)
        n_changed = (sub["DELTA_MAX_GSCORE"] != 0).sum()
        pct = 100.0 * n_changed / n if n else 0
        mean_d = sub.loc[sub["DELTA_MAX_GSCORE"] != 0, "DELTA_MAX_GSCORE"].mean() if n_changed else 0
        out[cls] = {"n": n, "pct": pct, "mean_delta": mean_d}
    return out


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    region_labels = []
    pct_from, pct_to = [], []
    mean_from, mean_to = [], []
    n_from, n_to = [], []

    for meta in REGIONS.values():
        s = load_direction_stats(meta["summary"])
        region_labels.append(meta["label"])
        pct_from.append(s["from_G"]["pct"])
        pct_to.append(s["to_G"]["pct"])
        mean_from.append(s["from_G"]["mean_delta"])
        mean_to.append(s["to_G"]["mean_delta"])
        n_from.append(s["from_G"]["n"])
        n_to.append(s["to_G"]["n"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    x = np.arange(len(region_labels))
    offset = 0.15

    # Left: % changed lollipop
    for i in range(len(x)):
        ax1.plot([x[i] - offset, x[i] - offset], [0, pct_from[i]], color=CLR_FROM, linewidth=2.2, zorder=2)
        ax1.plot([x[i] + offset, x[i] + offset], [0, pct_to[i]], color=CLR_TO, linewidth=2.2, zorder=2)
    ax1.scatter(x - offset, pct_from, color=CLR_FROM, s=100, zorder=3, label="from_G (G \u2192 non-G)")
    ax1.scatter(x + offset, pct_to, color=CLR_TO, s=100, zorder=3, label="to_G (non-G \u2192 G)")

    for i in range(len(x)):
        ax1.text(x[i] - offset, pct_from[i] + 2, f"n={n_from[i]}", ha="center", fontsize=8, color="#555")
        ax1.text(x[i] + offset, pct_to[i] + 2, f"n={n_to[i]}", ha="center", fontsize=8, color="#555")

    ax1.set_xticks(x, region_labels, fontsize=10)
    ax1.set_ylabel("% of variants with\nmax G-score change", fontsize=10)
    ax1.set_title("How often does G-score change?", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, max(max(pct_from), max(pct_to)) * 1.3)
    ax1.legend(fontsize=9, frameon=True)
    ax1.grid(axis="y", alpha=0.25)

    # Right: mean delta lollipop
    for i in range(len(x)):
        ax2.plot([x[i] - offset, x[i] - offset], [0, mean_from[i]], color=CLR_FROM, linewidth=2.2, zorder=2)
        ax2.plot([x[i] + offset, x[i] + offset], [0, mean_to[i]], color=CLR_TO, linewidth=2.2, zorder=2)
    ax2.scatter(x - offset, mean_from, color=CLR_FROM, s=100, zorder=3, label="from_G (G \u2192 non-G)")
    ax2.scatter(x + offset, mean_to, color=CLR_TO, s=100, zorder=3, label="to_G (non-G \u2192 G)")

    ax2.axhline(0, color="#333", linewidth=0.8)
    ax2.set_xticks(x, region_labels, fontsize=10)
    ax2.set_ylabel("Mean \u0394 max G-score\n(among changed variants)", fontsize=10)
    ax2.set_title("Direction and magnitude of change", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9, frameon=True)
    ax2.grid(axis="y", alpha=0.25)

    fig.suptitle("Mutation Direction and G-Quadruplex Disruption",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()

    outfile = OUTPUT_DIR / "combined_mutation_direction_lollipop.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
