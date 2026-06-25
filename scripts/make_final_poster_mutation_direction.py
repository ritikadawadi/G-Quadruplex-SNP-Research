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
        "label": "Upstream\nRegion",
        "summary": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
    },
    "3utr": {
        "label": "3\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
    },
}

CLR_FROM = "#d7301f"
CLR_TO = "#2171b5"


def load_direction_stats(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    rows = []
    for cls in ["from_G", "to_G"]:
        sub = df[df["G_CHANGE_CLASS"] == cls]
        n = len(sub)
        n_changed = (sub["DELTA_MAX_GSCORE"] != 0).sum()
        pct_changed = 100.0 * n_changed / n if n else 0
        mean_delta = sub.loc[sub["DELTA_MAX_GSCORE"] != 0, "DELTA_MAX_GSCORE"].mean() if n_changed else 0
        rows.append({
            "direction": cls,
            "n": n,
            "n_changed": n_changed,
            "pct_changed": pct_changed,
            "mean_delta": mean_delta,
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    region_keys = list(REGIONS.keys())
    for ax_idx, key in enumerate(region_keys):
        meta = REGIONS[key]
        stats = load_direction_stats(meta["summary"])
        from_g = stats[stats["direction"] == "from_G"].iloc[0]
        to_g = stats[stats["direction"] == "to_G"].iloc[0]

        labels = [f"from_G\n(G \u2192 non-G)\nn={int(from_g['n'])}",
                  f"to_G\n(non-G \u2192 G)\nn={int(to_g['n'])}"]

        # Left panel: % changed
        if ax_idx == 0:
            ax = axes[0]
            x = np.arange(3)
            width = 0.35
            # We'll build grouped bars per region
            break

    # Rebuild as grouped bar charts
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    region_labels = [REGIONS[k]["label"] for k in region_keys]
    x = np.arange(len(region_labels))
    width = 0.32

    pct_from, pct_to = [], []
    mean_from, mean_to = [], []
    n_from, n_to = [], []

    for key in region_keys:
        stats = load_direction_stats(REGIONS[key]["summary"])
        fg = stats[stats["direction"] == "from_G"].iloc[0]
        tg = stats[stats["direction"] == "to_G"].iloc[0]
        pct_from.append(fg["pct_changed"])
        pct_to.append(tg["pct_changed"])
        mean_from.append(fg["mean_delta"])
        mean_to.append(tg["mean_delta"])
        n_from.append(int(fg["n"]))
        n_to.append(int(tg["n"]))

    bars1 = ax1.bar(x - width/2, pct_from, width, color=CLR_FROM, label="from_G (G \u2192 non-G)", edgecolor="white")
    bars2 = ax1.bar(x + width/2, pct_to, width, color=CLR_TO, label="to_G (non-G \u2192 G)", edgecolor="white")

    for i in range(len(x)):
        ax1.text(x[i] - width/2, pct_from[i] + 1.5, f"n={n_from[i]}", ha="center", fontsize=8, color="#555")
        ax1.text(x[i] + width/2, pct_to[i] + 1.5, f"n={n_to[i]}", ha="center", fontsize=8, color="#555")

    ax1.set_xticks(x, region_labels, fontsize=10)
    ax1.set_ylabel("% of variants with\nmax G-score change", fontsize=10)
    ax1.set_title("How often does G-score change?", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, max(max(pct_from), max(pct_to)) * 1.25)
    ax1.legend(fontsize=9, frameon=True)
    ax1.grid(axis="y", alpha=0.25)

    bars3 = ax2.bar(x - width/2, mean_from, width, color=CLR_FROM, label="from_G (G \u2192 non-G)", edgecolor="white")
    bars4 = ax2.bar(x + width/2, mean_to, width, color=CLR_TO, label="to_G (non-G \u2192 G)", edgecolor="white")

    ax2.set_xticks(x, region_labels, fontsize=10)
    ax2.set_ylabel("Mean \u0394 max G-score\n(among changed variants)", fontsize=10)
    ax2.set_title("Direction and magnitude of change", fontsize=11, fontweight="bold")
    ax2.axhline(0, color="#333", linewidth=0.8)
    ax2.legend(fontsize=9, frameon=True)
    ax2.grid(axis="y", alpha=0.25)

    fig.suptitle("Mutation Direction and G-Quadruplex Disruption",
                 fontsize=13, fontweight="bold", y=1.01)

    fig.tight_layout()
    outfile = OUTPUT_DIR / "combined_mutation_direction.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
