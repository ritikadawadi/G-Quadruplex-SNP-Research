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
    "tss": {
        "label": "Upstream Region",
        "table": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "bio_clinical_cosmic_role_gscore_impact.tsv",
        "color": "#e6550d",
    },
    "5utr": {
        "label": "5'UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "bio_clinical_cosmic_role_gscore_impact.tsv",
        "color": "#2ca25f",
    },
    "3utr": {
        "label": "3'UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "bio_clinical_cosmic_role_gscore_impact.tsv",
        "color": "#756bb1",
    },
}

ROLES = ["TSG", "fusion", "oncogene"]


def load_region_df(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df[df["ROLE_SHORT"].isin(ROLES)].copy()


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    region_dfs = {key: load_region_df(meta["table"]).set_index("ROLE_SHORT") for key, meta in REGIONS.items()}

    x = np.arange(len(ROLES))
    n_regions = len(REGIONS)
    width = 0.25

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.45))
    fig.patch.set_facecolor("white")

    # Panel 1: percent changed
    ax = axes[0]
    for idx, (key, meta) in enumerate(REGIONS.items()):
        vals = [float(region_dfs[key].loc[role, "pct_gscore_changed"]) for role in ROLES]
        ns = [int(region_dfs[key].loc[role, "n_variants"]) for role in ROLES]
        offset = (idx - (n_regions - 1) / 2) * width
        xpos = x + offset
        bars = ax.bar(xpos, vals, width=width, color=meta["color"], label=meta["label"], edgecolor="white", linewidth=0.8)
        for bar, n in zip(bars, ns):
            bh = bar.get_height()
            if bh > 5:
                ax.text(bar.get_x() + bar.get_width() / 2, bh / 2, f"n={n}",
                        ha="center", va="center", fontsize=9, color="white", fontweight="bold")
            else:
                ax.text(bar.get_x() + bar.get_width() / 2, bh + 0.8, f"n={n}",
                        ha="center", va="bottom", fontsize=9, color="#333")

    ax.set_xticks(x, ["Tumor Suppressor", "Fusion Gene", "Oncogene"])
    ax.set_ylim(0, 75)
    ax.set_ylabel("% of variants with max G-score change", fontsize=16)
    ax.set_title(
        "How often does the strongest\npredicted G-score change?",
        fontsize=20,
        fontweight="bold",
    )
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=14)
    ax.legend(frameon=True, loc="upper right", fontsize=13)
    ax.grid(axis="y", alpha=0.25)
    ax.text(0.97, 0.78, "n = no. of 'G'-related\nSNVs per category",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color="#888888",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))

    # Panel 2: magnitude among changed
    region_dfs = {key: load_region_df(meta["table"]).set_index("ROLE_SHORT") for key, meta in REGIONS.items()}

    ax = axes[1]
    for idx, (key, meta) in enumerate(REGIONS.items()):
        vals = [float(region_dfs[key].loc[role, "mean_abs_delta_when_changed"]) for role in ROLES]
        offset = (idx - (n_regions - 1) / 2) * width
        xpos = x + offset
        ax.bar(xpos, vals, width=width, color=meta["color"], label=meta["label"], edgecolor="white", linewidth=0.8)

    ax.set_xticks(x, ["Tumor Suppressor", "Fusion Gene", "Oncogene"])
    ax.set_ylabel("Mean |Δ max G-score| when changed", fontsize=16)
    ax.set_title("When it changes, how large is the shift?", fontsize=20, fontweight="bold")
    ax.set_ylim(0, 14)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=14)
    ax.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Cancer gene role and G-score disruption across regions",
        fontsize=22,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(
        left=0.07,
        right=0.99,
        bottom=0.10,
        top=0.78,
        wspace=0.30,
    )

    outfile = OUTPUT_DIR / "combined_cosmic_role_g4_change_5utr_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
