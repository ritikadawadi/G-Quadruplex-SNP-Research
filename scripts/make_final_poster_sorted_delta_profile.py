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


PROJECT_ROOT = Path(__file__).resolve().parent.parent
COMPARE_ROOT = PROJECT_ROOT / "compare_mutant_wild"
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"


REGIONS = {
    "tss": {
        "label": "Upstream Region",
        "summary": COMPARE_ROOT / "tss" / "tss_upstream500.comparison_summary.tsv",
        "color": "#e6550d",
    },
    "5utr": {
        "label": "5'UTR",
        "summary": COMPARE_ROOT / "5utr" / "5utr.comparison_summary.tsv",
        "color": "#08519c",
    },
    "3utr": {
        "label": "3'UTR",
        "summary": COMPARE_ROOT / "3utr" / "3utr.comparison_summary.tsv",
        "color": "#31a354",
    },
}


def load_sorted_delta(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    return pd.to_numeric(df["DELTA_MAX_GSCORE"], errors="coerce").fillna(0).sort_values().reset_index(drop=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    fig, ax = plt.subplots(figsize=(11, 6))
    for meta in REGIONS.values():
        values = load_sorted_delta(meta["summary"])
        ax.plot(
            range(1, len(values) + 1),
            values.values,
            color=meta["color"],
            linewidth=2.2,
            label=f"{meta['label']} (n={len(values)})",
        )

    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(
        "Mutation-Caused Changes in the G-Score in UTR and Upstream Regions of Cancer Genes",
        fontweight="bold",
        fontsize=20,
        pad=12,
    )
    ax.set_xlabel("Mutant SNPs", fontsize=16)
    ax.set_ylabel("Delta max G-score", fontsize=16)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(frameon=True, fontsize=14)
    ax.text(
        0.97, 0.03,
        "n = no. of mutant SNP\nsequences per region",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=14, color="#555555",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9),
    )
    fig.tight_layout()

    outfile = OUTPUT_DIR / "combined_11_sorted_delta_max_gscore_profile_5utr_tss_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
