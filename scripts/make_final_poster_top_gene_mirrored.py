#!/usr/bin/env python3
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

REGIONS = {
    "5utr": {
        "label": "5\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
    },
    "tss": {
        "label": "TSS",
        "summary": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
    },
    "3utr": {
        "label": "3\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
    },
}

N_LOSS = 3
N_GAIN = 2

CLR_LOSS = "#d7301f"
CLR_GAIN = "#2171b5"


def load_mixed(summary_tsv: Path) -> pd.DataFrame:
    """Top losses + top gains from the full comparison summary."""
    df = pd.read_csv(summary_tsv, sep="\t")
    g = df.groupby("GENE").agg(
        value=("DELTA_MAX_GSCORE", "sum"),
        n=("DELTA_MAX_GSCORE", "count"),
    ).reset_index()
    g["abs_val"] = g["value"].abs()
    losses = g[g["value"] < 0].sort_values("abs_val", ascending=False).head(N_LOSS)
    gains = g[g["value"] > 0].sort_values("value", ascending=False).head(N_GAIN)
    combined = pd.concat([losses, gains]).sort_values("value")
    return combined[["GENE", "value", "n"]]


def panel(ax, df, *, title):
    work = df.copy().reset_index(drop=True)
    y = np.arange(len(work))
    vals = work["value"].astype(float).values
    colors = [CLR_LOSS if v < 0 else CLR_GAIN for v in vals]

    ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.8, height=0.65)
    ax.axvline(0, color="#333333", linewidth=0.9)

    mag = max(abs(vals.min()), abs(vals.max())) if len(vals) else 1.0
    xlim = mag * 1.30
    ax.set_xlim(-xlim, xlim)
    ax.set_yticks(y, work["GENE"].tolist())
    ax.invert_yaxis()
    ax.set_title(title, fontsize=11, pad=4)
    ax.set_xlabel("Net \u0394 max G-score (sum across SNVs)", fontsize=9)
    ax.tick_params(labelsize=9)
    ax.grid(axis="x", alpha=0.25)

    for i, row in work.iterrows():
        val = float(row["value"])
        n_snvs = int(row["n"])
        offset = xlim * 0.03
        if val < 0:
            tx, ha = val - offset, "right"
        else:
            tx, ha = val + offset, "left"
        ax.text(tx, y[i], f"n={n_snvs}", ha=ha, va="center",
                fontsize=8, color="#555555")

    ax.text(0.97, 0.97, "n = no. of 'G'-related\nSNVs per gene",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7, color="#888888",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    signed_5 = load_mixed(REGIONS["5utr"]["summary"])
    signed_tss = load_mixed(REGIONS["tss"]["summary"])
    signed_3 = load_mixed(REGIONS["3utr"]["summary"])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.patch.set_facecolor("white")

    panel(axes[0], signed_5, title="5\u2032UTR")
    panel(axes[1], signed_tss, title="Upstream Region")
    panel(axes[2], signed_3, title="3\u2032UTR")

    fig.suptitle(
        "Top genes with strongest net G-score changes",
        fontsize=13, fontweight="bold", y=1.02,
    )

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=CLR_LOSS, label="G4 loss (negative \u0394)"),
        Patch(facecolor=CLR_GAIN, label="G4 gain (positive \u0394)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center",
               ncol=2, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, -0.04))

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    outfile = OUTPUT_DIR / "combined_top_gene_mirrored_5utr_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
