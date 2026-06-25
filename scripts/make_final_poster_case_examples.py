#!/usr/bin/env python3
"""Case-study figure: one 5'UTR and one 3'UTR example of G4 disruption."""
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"


CASES = [
    {
        "region": "5\u2032UTR",
        "gene": "MSH2",
        "role": "TSG",
        "variant": "chr2:47,403,183  G \u2192 C",
        "clinvar": "ClinVar 182605",
        "disease": "Lynch syndrome / Hereditary nonpolyposis colorectal cancer",
        "wt_seq": "GGAGGTGAGGAGG",
        "mut_seq": "\u2014  (motif destroyed)",
        "wt_gscore": 34,
        "mut_gscore": 0,
        "delta": -34,
        "wt_qgrs": 1,
        "mut_qgrs": 0,
        "effect": "Complete G4 motif loss",
        "accent": "#2ca25f",
        "bg": "#f0faf4",
    },
    {
        "region": "3\u2032UTR",
        "gene": "SMAD3",
        "role": "TSG",
        "variant": "chr15:67,191,109  G \u2192 T",
        "clinvar": "ClinVar 316878",
        "disease": "Loeys-Dietz syndrome / Familial thoracic aortic aneurysm",
        "wt_seq": "GGGATGGGGAGGGAGGG",
        "mut_seq": "GGTGATGGTGGGGAGATGATGGGCTAAACAGG",
        "wt_gscore": 63,
        "mut_gscore": 36,
        "delta": -27,
        "wt_qgrs": 21,
        "mut_qgrs": 21,
        "effect": "G4 score weakened (motif retained)",
        "accent": "#756bb1",
        "bg": "#f5f0fa",
    },
]


def draw_case(ax, case, y_start):
    """Draw one case study card."""
    accent = case["accent"]
    bg = case["bg"]

    card_w, card_h = 0.92, 0.40
    card_x = 0.04
    card_y = y_start

    patch = FancyBboxPatch(
        (card_x, card_y), card_w, card_h,
        boxstyle="round,pad=0.01,rounding_size=0.015",
        linewidth=1.8, edgecolor=accent, facecolor=bg,
    )
    ax.add_patch(patch)

    # Region + Gene header
    hdr_y = card_y + card_h - 0.04
    ax.text(card_x + 0.03, hdr_y, f"{case['region']}  |  ",
            fontsize=13, fontweight="bold", color="#333333", va="top",
            fontfamily="sans-serif")
    ax.text(card_x + 0.19, hdr_y, f"{case['gene']}",
            fontsize=14, fontweight="bold", color=accent, va="top")
    ax.text(card_x + 0.27, hdr_y, f"({case['role']})",
            fontsize=11, color="#666666", va="top")

    ax.text(card_x + 0.50, hdr_y, case["variant"],
            fontsize=10, color="#333333", va="top",
            fontfamily="monospace")
    ax.text(card_x + 0.78, hdr_y, case["clinvar"],
            fontsize=9, color="#888888", va="top")

    # Disease
    dis_y = hdr_y - 0.05
    ax.text(card_x + 0.03, dis_y, "Disease:",
            fontsize=9, fontweight="bold", color="#555555", va="top")
    ax.text(card_x + 0.12, dis_y, case["disease"],
            fontsize=9, color="#555555", va="top", style="italic")

    # Separator
    ax.plot([card_x + 0.02, card_x + card_w - 0.02],
            [dis_y - 0.035, dis_y - 0.035],
            color=accent, lw=0.6, alpha=0.4)

    # Wild-type row
    row1_y = dis_y - 0.065
    ax.text(card_x + 0.03, row1_y, "Wild-type",
            fontsize=10, fontweight="bold", color="#1a6b3c", va="top")
    ax.text(card_x + 0.16, row1_y, f"G-score: {case['wt_gscore']}",
            fontsize=10, color="#333333", va="top")
    ax.text(card_x + 0.30, row1_y, f"QGRS count: {case['wt_qgrs']}",
            fontsize=10, color="#333333", va="top")
    ax.text(card_x + 0.48, row1_y, f"Top motif:  {case['wt_seq']}",
            fontsize=8.5, color="#333333", va="top",
            fontfamily="monospace")

    # Mutant row
    row2_y = row1_y - 0.05
    ax.text(card_x + 0.03, row2_y, "Mutant",
            fontsize=10, fontweight="bold", color="#b33a2a", va="top")
    ax.text(card_x + 0.16, row2_y, f"G-score: {case['mut_gscore']}",
            fontsize=10, color="#333333", va="top")
    ax.text(card_x + 0.30, row2_y, f"QGRS count: {case['mut_qgrs']}",
            fontsize=10, color="#333333", va="top")
    mut_seq_text = case["mut_seq"]
    if len(mut_seq_text) > 35:
        mut_seq_text = mut_seq_text[:34] + "…"
    ax.text(card_x + 0.48, row2_y, f"Top motif:  {mut_seq_text}",
            fontsize=8.5, color="#333333", va="top",
            fontfamily="monospace")

    # Delta summary row
    row3_y = row2_y - 0.055
    ax.plot([card_x + 0.02, card_x + card_w - 0.02],
            [row3_y + 0.02, row3_y + 0.02],
            color=accent, lw=0.6, alpha=0.4)

    delta_color = "#c0392b" if case["delta"] < 0 else "#27ae60"
    ax.text(card_x + 0.03, row3_y,
            f"\u0394 G-score: {case['delta']:+d}",
            fontsize=11, fontweight="bold", color=delta_color, va="top")
    ax.text(card_x + 0.25, row3_y,
            f"\u0394 QGRS count: {case['mut_qgrs'] - case['wt_qgrs']:+d}",
            fontsize=11, fontweight="bold", color=delta_color if case['mut_qgrs'] != case['wt_qgrs'] else "#555555", va="top")
    ax.text(card_x + 0.52, row3_y,
            case["effect"],
            fontsize=10, fontweight="bold", color=accent, va="top")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(0.5, 0.96,
            "Case examples: G-related SNVs disrupting predicted G4 motifs",
            ha="center", va="top", fontsize=15, fontweight="bold",
            color="#1e293b")

    draw_case(ax, CASES[0], y_start=0.52)
    draw_case(ax, CASES[1], y_start=0.06)

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = OUTPUT_DIR / "case_examples_5utr_3utr.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
