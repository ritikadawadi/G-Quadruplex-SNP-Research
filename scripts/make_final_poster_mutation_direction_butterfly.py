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
        out[cls] = {"n": n, "n_changed": n_changed, "pct": pct, "mean_delta": mean_d}
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    y = np.arange(len(region_labels))
    height = 0.55

    # --- Panel 1: % changed butterfly ---
    ax1.barh(y, [-p for p in pct_from], height, color=CLR_FROM, edgecolor="white",
             linewidth=0.8, label="from_G (G \u2192 non-G)")
    ax1.barh(y, pct_to, height, color=CLR_TO, edgecolor="white",
             linewidth=0.8, label="to_G (non-G \u2192 G)")
    ax1.axvline(0, color="#333", linewidth=1)

    max_pct = max(max(pct_from), max(pct_to))
    ax1.set_xlim(-max_pct * 1.45, max_pct * 1.45)

    for i in range(len(y)):
        ax1.text(-pct_from[i] - max_pct * 0.04, y[i], f"{pct_from[i]:.1f}%\n(n={n_from[i]})",
                 ha="right", va="center", fontsize=8.5, color=CLR_FROM, fontweight="bold")
        ax1.text(pct_to[i] + max_pct * 0.04, y[i], f"{pct_to[i]:.1f}%\n(n={n_to[i]})",
                 ha="left", va="center", fontsize=8.5, color=CLR_TO, fontweight="bold")

    ax1.set_yticks(y, region_labels, fontsize=11)
    ax1.set_xlabel("% of variants with max G-score change", fontsize=10)
    ax1.set_title("How often does G-score change?", fontsize=11, fontweight="bold")
    ax1.legend(loc="lower left", fontsize=8.5, frameon=True)
    ax1.grid(axis="x", alpha=0.25)
    ax1.invert_yaxis()
    ticks = ax1.get_xticks()
    ax1.set_xticklabels([f"{abs(t):.0f}%" for t in ticks])

    # --- Panel 2: mean delta butterfly ---
    ax2.barh(y, mean_from, height, color=CLR_FROM, edgecolor="white",
             linewidth=0.8, label="from_G (G \u2192 non-G)")
    ax2.barh(y, mean_to, height, color=CLR_TO, edgecolor="white",
             linewidth=0.8, label="to_G (non-G \u2192 G)")
    ax2.axvline(0, color="#333", linewidth=1)

    max_abs = max(max(abs(m) for m in mean_from), max(abs(m) for m in mean_to))
    ax2.set_xlim(-max_abs * 1.35, max_abs * 1.35)

    for i in range(len(y)):
        if mean_from[i] < 0:
            ax2.text(mean_from[i] - max_abs * 0.04, y[i], f"{mean_from[i]:.1f}",
                     ha="right", va="center", fontsize=9, color=CLR_FROM, fontweight="bold")
        else:
            ax2.text(mean_from[i] + max_abs * 0.04, y[i], f"{mean_from[i]:.1f}",
                     ha="left", va="center", fontsize=9, color=CLR_FROM, fontweight="bold")
        ax2.text(mean_to[i] + max_abs * 0.04, y[i], f"+{mean_to[i]:.1f}",
                 ha="left", va="center", fontsize=9, color=CLR_TO, fontweight="bold")

    ax2.set_yticks(y, region_labels, fontsize=11)
    ax2.set_xlabel("Mean \u0394 max G-score (among changed variants)", fontsize=10)
    ax2.set_title("Direction and magnitude of change", fontsize=11, fontweight="bold")
    ax2.legend(loc="lower left", fontsize=8.5, frameon=True)
    ax2.grid(axis="x", alpha=0.25)
    ax2.invert_yaxis()

    fig.suptitle("Mutation Direction and G-Quadruplex Disruption",
                 fontsize=14, fontweight="bold", y=1.02)

    fig.text(0.5, -0.03,
             "from_G: the reference allele is G and mutation removes it (G4 loss expected). "
             "to_G: the mutant allele introduces a G (G4 gain expected).",
             ha="center", fontsize=9, color="#666")

    fig.tight_layout()
    outfile = OUTPUT_DIR / "combined_mutation_direction_butterfly.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
