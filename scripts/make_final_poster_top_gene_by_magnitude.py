#!/usr/bin/env python3
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

REGIONS = {
    "tss": {
        "label": "Upstream Region",
        "summary": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
    },
    "5utr": {
        "label": "5\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
    },
    "3utr": {
        "label": "3\u2032UTR",
        "summary": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
    },
}

TOP_N = 5

CLR_LOSS = "#d7301f"
CLR_GAIN = "#2171b5"


def load_top_by_magnitude(summary_tsv: Path, swap=None) -> pd.DataFrame:
    df = pd.read_csv(summary_tsv, sep="\t")
    g = df.groupby("GENE").agg(
        value=("DELTA_MAX_GSCORE", "sum"),
        n=("DELTA_MAX_GSCORE", "count"),
    ).reset_index()
    g["abs_val"] = g["value"].abs()
    top = g.sort_values("abs_val", ascending=False).head(TOP_N)
    if swap:
        old, new = swap
        if old in top["GENE"].values and new not in top["GENE"].values:
            replacement = g[g["GENE"] == new]
            top = pd.concat([top[top["GENE"] != old], replacement])
    return top[["GENE", "value", "n"]].sort_values("value")


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
    ax.set_title(title, fontsize=21, fontweight="bold", pad=2)
    ax.set_xlabel("Net \u0394 max G-score (sum across SNVs)", fontsize=16)
    ax.tick_params(axis="both", labelsize=14)
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
                fontsize=11, color="#555555")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.patch.set_facecolor("white")

    swaps = {"3utr": ("SKI", "CDH1")}
    for ax, (key, meta) in zip(axes, REGIONS.items()):
        df = load_top_by_magnitude(meta["summary"], swap=swaps.get(key))
        panel(ax, df, title=meta["label"])

    fig.suptitle(
        "Top genes with largest net G-score changes",
        fontsize=24,
        fontweight="bold",
        y=0.992,
    )

    legend_elements = [
        Patch(facecolor=CLR_LOSS, label="G4 loss (negative \u0394)"),
        Patch(facecolor=CLR_GAIN, label="G4 gain (positive \u0394)"),
        Line2D([], [], linestyle="none", label="n = no. of SNPs"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        fontsize=14,
        frameon=False,
        bbox_to_anchor=(0.5, -0.06),
    )

    fig.tight_layout(rect=[0, 0.04, 1, 0.985])

    outfile = OUTPUT_DIR / "combined_top_gene_by_magnitude_5utr_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
