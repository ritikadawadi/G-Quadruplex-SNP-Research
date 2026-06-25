#!/usr/bin/env python3
"""Generate a poster-ready summary statistics table image."""
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results_compare_figures" / "final poster figures"

HEADER_BG = "#6D9BC3"
HEADER_FG = "white"
SECTION_BG = "#D6E6F0"
ROW_BG_ALT = ("#FFFFFF", "#F4F8FB")
INDENT = "    "
# Table width and column split (data coords): wider layout + clear two-column gutter
TABLE_W = 9.0
COL_SPLIT = 6.35
LABEL_X = 0.35
VALUE_X = TABLE_W - 0.35


def build_table():
    rows = [
        # (label, value, is_section_header)
        ("Pipeline Input", "", True),
        ("COSMIC Cancer Census Genes", "763", False),
        ("Upstream Region (500 bp)", "", True),
        (f"{INDENT}All ClinVar–QGRS Overlaps", "442", False),
        (f"{INDENT}Variants with Δ max G-score ≠ 0", "62  (20.2%)", False),
        (f"{INDENT}Wild-type / Mutant QGRS Motifs", "1,048 / 1,023", False),
        ("5′UTR", "", True),
        (f"{INDENT}All ClinVar–QGRS Overlaps", "602", False),
        (f"{INDENT}Variants with Δ max G-score ≠ 0", "144  (35.5%)", False),
        (f"{INDENT}Wild-type / Mutant QGRS Motifs", "832 / 786", False),
        ("3′UTR", "", True),
        (f"{INDENT}All ClinVar–QGRS Overlaps", "974", False),
        (f"{INDENT}Variants with Δ max G-score ≠ 0", "54  (8.5%)", False),
        (f"{INDENT}Wild-type / Mutant QGRS Motifs", "7,845 / 7,760", False),
    ]
    return rows


def render_table(rows, outfile):
    n = len(rows)
    fig_h = 0.52 * n + 0.65
    fig_w = 6 * (TABLE_W / 7.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, TABLE_W)
    ax.set_ylim(0, n)
    ax.axis("off")

    alt_idx = 0
    for i, (label, value, is_header) in enumerate(reversed(rows)):
        y = i
        if not label and is_header:
            continue
        if is_header:
            ax.add_patch(plt.Rectangle((0, y), TABLE_W, 1, fc=SECTION_BG, ec="none"))
            ax.text(LABEL_X, y + 0.5, label, fontsize=18, fontweight="bold",
                    va="center", ha="left", color="#1a1a1a")
        else:
            bg = ROW_BG_ALT[alt_idx % 2]
            alt_idx += 1
            ax.add_patch(plt.Rectangle((0, y), TABLE_W, 1, fc=bg, ec="none"))
            ax.text(LABEL_X, y + 0.5, label, fontsize=17, va="center", ha="left", color="#222")
            ax.text(VALUE_X, y + 0.5, value, fontsize=17, va="center", ha="right",
                    fontweight="bold", color="#222")

    for i in range(n + 1):
        ax.plot([0, TABLE_W], [i, i], color="#bbb", linewidth=0.5)
    ax.plot([0, 0], [0, n], color="#bbb", linewidth=0.5)
    ax.plot([TABLE_W, TABLE_W], [0, n], color="#bbb", linewidth=0.5)
    ax.plot(
        [COL_SPLIT, COL_SPLIT],
        [0, n],
        color="#c8c8c8",
        linewidth=0.75,
        linestyle="-",
    )

    fig.tight_layout(pad=0.3)
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = build_table()
    render_table(rows, OUTPUT_DIR / "project_summary_table.png")


if __name__ == "__main__":
    main()
