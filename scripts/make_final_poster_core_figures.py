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
import textwrap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

REGIONS = {
    "tss": {
        "label": "TSS",
        "path": PROJECT_ROOT / "compare_mutant_wild" / "tss" / "tss_upstream500.comparison_summary.tsv",
        "color": "#3182bd",
    },
    "5utr": {
        "label": "5'UTR",
        "path": PROJECT_ROOT / "compare_mutant_wild" / "5utr" / "5utr.comparison_summary.tsv",
        "color": "#31a354",
    },
    "3utr": {
        "label": "3'UTR",
        "path": PROJECT_ROOT / "compare_mutant_wild" / "3utr" / "3utr.comparison_summary.tsv",
        "color": "#756bb1",
    },
}

PIPELINE_COUNTS = {
    "genes": 763,
    "tss_sequences": 756,
    "utr5_sequences": 753,
    "utr3_sequences": 753,
    "all_overlaps": 2018,
    "g_related": 1351,
    "unique_genes": 259,
}


def load_region_stats() -> pd.DataFrame:
    rows = []
    for key, meta in REGIONS.items():
        df = pd.read_csv(meta["path"], sep="\t")
        delta_g = pd.to_numeric(df["DELTA_MAX_GSCORE"], errors="coerce").fillna(0)
        delta_q = pd.to_numeric(df["DELTA_QGRS_COUNT"], errors="coerce").fillna(0)
        rows.append(
            {
                "region": meta["label"],
                "color": meta["color"],
                "n_variants": len(df),
                "pct_gscore_changed": 100 * (delta_g != 0).mean(),
                "mean_abs_delta_when_changed": delta_g[delta_g != 0].abs().mean() if (delta_g != 0).any() else 0.0,
                "pct_top_sequence_changed": 100 * df["TOP_SEQUENCE_CHANGED"].astype(str).eq("y").mean(),
                "pct_gscore_lower": 100 * (delta_g < 0).mean(),
                "pct_gscore_same": 100 * (delta_g == 0).mean(),
                "pct_gscore_higher": 100 * (delta_g > 0).mean(),
                "pct_qgrs_fewer": 100 * (delta_q < 0).mean(),
                "pct_qgrs_same": 100 * (delta_q == 0).mean(),
                "pct_qgrs_more": 100 * (delta_q > 0).mean(),
            }
        )
    return pd.DataFrame(rows)


def draw_box(ax, xy, width, height, title, subtitle, facecolor):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.2,
        edgecolor="#3f4a56",
        facecolor=facecolor,
    )
    ax.add_patch(box)
    x, y = xy
    title_text = textwrap.fill(title, width=18)
    subtitle_text = textwrap.fill(subtitle, width=22)
    ax.text(
        x + width / 2,
        y + height * 0.63,
        title_text,
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color="#1f1f1f",
        linespacing=1.1,
        clip_on=True,
    )
    ax.text(
        x + width / 2,
        y + height * 0.28,
        subtitle_text,
        ha="center",
        va="center",
        fontsize=8.8,
        color="#374151",
        linespacing=1.1,
        clip_on=True,
    )


def make_pipeline_overview():
    fig, ax = plt.subplots(figsize=(14, 5.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle(
        "Project workflow: from cancer genes to predicted G4-disrupting variants",
        fontsize=16,
        fontweight="bold",
        y=0.94,
    )

    boxes = [
        ((0.02, 0.34), 0.13, 0.28, "COSMIC cancer genes", f"{PIPELINE_COUNTS['genes']} genes", "#dbeafe"),
        ((0.22, 0.34), 0.13, 0.28, "Representative transcripts", "1 transcript per gene", "#e0f2fe"),
        ((0.42, 0.34), 0.15, 0.28, "Regulatory regions", "TSS, 5'UTR, 3'UTR\nstrand-aware extraction", "#ecfccb"),
        ((0.64, 0.34), 0.13, 0.28, "ClinVar overlap", f"{PIPELINE_COUNTS['all_overlaps']} overlaps", "#fef3c7"),
        ((0.84, 0.34), 0.13, 0.28, "Wild vs mutant QGRS", f"{PIPELINE_COUNTS['g_related']} G-related SNVs\n{PIPELINE_COUNTS['unique_genes']} genes retained", "#fee2e2"),
    ]
    for xy, w, h, title, subtitle, fc in boxes:
        draw_box(ax, xy, w, h, title, subtitle, fc)

    for i in range(len(boxes) - 1):
        x1 = boxes[i][0][0] + boxes[i][1]
        x2 = boxes[i + 1][0][0]
        y = 0.48
        arrow = FancyArrowPatch((x1 + 0.004, y), (x2 - 0.004, y), arrowstyle="-|>", mutation_scale=20, linewidth=1.8, color="#4b5563")
        ax.add_patch(arrow)

    ax.text(
        0.5,
        0.16,
        "Reference QGRS scans used 756 TSS, 753 5'UTR, and 753 3'UTR sequences. ClinVar variants were intersected with predicted QGRS loci, then filtered to transcript-oriented G-related SNVs before mutant-sequence rescanning.",
        ha="center",
        va="center",
        fontsize=10,
        color="#374151",
        wrap=True,
    )
    fig.tight_layout(rect=[0.01, 0.03, 0.99, 0.90])
    fig.savefig(OUTPUT_DIR / "pipeline_overview_workflow.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_cross_region_summary(stats: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.2))
    colors = stats["color"].tolist()
    x = np.arange(len(stats))
    labels = stats["region"].tolist()

    panels = [
        ("n_variants", "G-related SNVs compared", "Number of variant comparisons"),
        ("pct_gscore_changed", "Variants with max G-score change", "% with Δ max G-score ≠ 0"),
        ("mean_abs_delta_when_changed", "Magnitude when change happens", "Mean |Δ max G-score| among changers"),
        ("pct_top_sequence_changed", "Top motif sequence changed", "% with top-sequence change"),
    ]

    for ax, (col, title, ylabel) in zip(axes.flat, panels):
        vals = stats[col].astype(float).values
        bars = ax.bar(x, vals, color=colors, edgecolor="white", linewidth=0.8, width=0.62)
        ax.set_xticks(x, labels)
        ax.set_title(title, fontsize=12, pad=6)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.25)
        for bar, val in zip(bars, vals):
            fmt = f"{val:.1f}" if abs(val - round(val)) > 1e-6 else f"{int(round(val))}"
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(vals) * 0.03, fmt, ha="center", va="bottom", fontsize=9, color="#333333")

    fig.suptitle("Cross-region summary of clinically observed G4 effects", fontsize=16, fontweight="bold", y=0.96)
    fig.tight_layout(rect=[0.03, 0.04, 0.99, 0.93])
    fig.savefig(OUTPUT_DIR / "cross_region_summary.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_effect_direction_summary(stats: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.8))
    region_colors = stats["color"].tolist()
    y = np.arange(len(stats))
    loss_color = "#E68613"
    neutral_color = "#CFCFCF"
    gain_color = "#4472C4"

    # Max g-score direction
    ax = axes[0]
    lower = -stats["pct_gscore_lower"].astype(float).values
    same = stats["pct_gscore_same"].astype(float).values
    higher = stats["pct_gscore_higher"].astype(float).values
    ax.barh(y, lower, color=loss_color, edgecolor="white", linewidth=0.8, label="Lower max G-score")
    ax.barh(y, same, left=0, color=neutral_color, edgecolor="white", linewidth=0.8, label="Unchanged")
    ax.barh(y, higher, left=same, color=gain_color, edgecolor="white", linewidth=0.8, label="Higher max G-score")
    ax.axvline(0, color="#222222", linewidth=1)
    ax.set_yticks(y, stats["region"].tolist())
    ax.set_title("Direction of max G-score change", fontsize=12, pad=6)
    ax.set_xlabel("Percent of variants")
    ax.grid(axis="x", alpha=0.25)

    # QGRS count direction
    ax = axes[1]
    fewer = -stats["pct_qgrs_fewer"].astype(float).values
    same_q = stats["pct_qgrs_same"].astype(float).values
    more = stats["pct_qgrs_more"].astype(float).values
    ax.barh(y, fewer, color=loss_color, edgecolor="white", linewidth=0.8, label="Fewer QGRS motifs")
    ax.barh(y, same_q, left=0, color=neutral_color, edgecolor="white", linewidth=0.8, label="Unchanged")
    ax.barh(y, more, left=same_q, color=gain_color, edgecolor="white", linewidth=0.8, label="More QGRS motifs")
    ax.axvline(0, color="#222222", linewidth=1)
    ax.set_yticks(y, stats["region"].tolist())
    ax.set_title("Direction of QGRS count change", fontsize=12, pad=6)
    ax.set_xlabel("Percent of variants")
    ax.grid(axis="x", alpha=0.25)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.93))
    fig.suptitle("How do mutations shift predicted G4 outcomes?", fontsize=16, fontweight="bold", y=0.98)
    fig.text(0.5, 0.03, "Left-going orange bars show losses or decreases; right-going blue bars show gains or increases. Gray indicates no change.", ha="center", fontsize=10, color="#444444")
    fig.tight_layout(rect=[0.02, 0.07, 0.99, 0.88])
    fig.savefig(OUTPUT_DIR / "effect_direction_summary.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    stats = load_region_stats()
    make_pipeline_overview()
    make_cross_region_summary(stats)
    make_effect_direction_summary(stats)
    print("Wrote core final poster figures")


if __name__ == "__main__":
    main()
