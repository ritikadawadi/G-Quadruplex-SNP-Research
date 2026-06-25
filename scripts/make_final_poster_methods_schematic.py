#!/usr/bin/env python3
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


def rounded_box(ax, x, y, w, h, fc="#f0f4f8", ec="#cbd5e1", lw=1.2):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.01,rounding_size=0.015",
        linewidth=lw, edgecolor=ec, facecolor=fc,
    )
    ax.add_patch(patch)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.5, 9.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── Title ──
    ax.text(0.5, 0.96, "Methods", ha="center", va="top",
            fontsize=22, fontweight="bold", color="#1e293b")

    # ── Data sources row ──
    ax.text(0.5, 0.90, "Data Sources", ha="center", va="top",
            fontsize=13, fontweight="bold", color="#475569")

    sources = [
        {"cx": 0.13, "label_top": "COSMIC", "label_bot": "Cancer Gene\nCensus",
         "fc": "#dbeafe"},
        {"cx": 0.37, "label_top": "GENCODE +", "label_bot": "GRCh38",
         "fc": "#dcfce7"},
        {"cx": 0.63, "label_top": "QGRS", "label_bot": "Mapper",
         "fc": "#fef9c3"},
        {"cx": 0.87, "label_top": "ClinVar", "label_bot": "NCBI",
         "fc": "#fce7f3"},
    ]

    box_w, box_h = 0.18, 0.10
    box_y = 0.76
    for s in sources:
        bx = s["cx"] - box_w / 2
        rounded_box(ax, bx, box_y, box_w, box_h, fc=s["fc"], ec="#94a3b8")
        ax.text(s["cx"], box_y + box_h / 2 + 0.005, "(logo)",
                ha="center", va="center", fontsize=9, color="#9ca3af",
                style="italic")
        ax.text(s["cx"], box_y - 0.025, s["label_top"],
                ha="center", va="top", fontsize=10, fontweight="bold",
                color="#334155")
        ax.text(s["cx"], box_y - 0.045, s["label_bot"],
                ha="center", va="top", fontsize=9, color="#64748b")

    # ── Separator ──
    ax.plot([0.08, 0.92], [0.68, 0.68], color="#cbd5e1", lw=0.8)

    # ── Pipeline steps (numbered list) ──
    ax.text(0.5, 0.65, "Pipeline", ha="center", va="top",
            fontsize=13, fontweight="bold", color="#475569")

    steps = [
        ("I.", "Select one representative transcript per gene via MANE Select",
         "763 COSMIC cancer genes"),
        ("II.", "Extract regulatory regions from GRCh38",
         "TSS (500 bp upstream), 5\u2032UTR, 3\u2032UTR"),
        ("III.", "Scan reference sequences with QGRS Mapper",
         "Predict G4 motifs in wild-type"),
        ("IV.", "Overlap with ClinVar and filter to G-related SNVs",
         "n = 1,351 variants across all regions"),
        ("V.", "Generate mutant sequences (one SNV per FASTA)",
         "Substitute variant allele into reference"),
        ("VI.", "Re-scan mutant sequences with QGRS Mapper",
         "Compare wild-type vs mutant G4 predictions"),
    ]

    y_start = 0.60
    y_step = 0.065
    for i, (num, desc, detail) in enumerate(steps):
        y = y_start - i * y_step
        ax.text(0.10, y, num, ha="right", va="top",
                fontsize=11, fontweight="bold", color="#1e40af")
        ax.text(0.12, y, desc, ha="left", va="top",
                fontsize=10.5, fontweight="bold", color="#1e293b")
        ax.text(0.12, y - 0.023, detail, ha="left", va="top",
                fontsize=9, color="#64748b")

    # ── Separator ──
    sep_y = y_start - len(steps) * y_step + 0.015
    ax.plot([0.08, 0.92], [sep_y, sep_y], color="#cbd5e1", lw=0.8)

    # ── Analysis & Output ──
    ax.text(0.5, sep_y - 0.02, "Analysis & Output", ha="center", va="top",
            fontsize=13, fontweight="bold", color="#475569")

    outputs = [
        ("VII.", "Quantify G4 change",
         "\u0394 max G-score, \u0394 QGRS count, motif gain/loss"),
        ("VIII.", "Biological interpretation",
         "ClinVar significance, disease terms (CLNDN), COSMIC cancer roles"),
    ]

    y_out = sep_y - 0.06
    for i, (num, desc, detail) in enumerate(outputs):
        y = y_out - i * y_step
        ax.text(0.10, y, num, ha="right", va="top",
                fontsize=11, fontweight="bold", color="#1e40af")
        ax.text(0.12, y, desc, ha="left", va="top",
                fontsize=10.5, fontweight="bold", color="#1e293b")
        ax.text(0.12, y - 0.023, detail, ha="left", va="top",
                fontsize=9, color="#64748b")

    # ── Bottom summary ──
    bot_y = y_out - len(outputs) * y_step - 0.005
    ax.plot([0.08, 0.92], [bot_y, bot_y], color="#cbd5e1", lw=0.8)
    ax.text(0.5, bot_y - 0.02,
            "Output: figures and tables identifying which clinically observed\n"
            "G-related variants alter predicted regulatory G4 motifs in cancer genes",
            ha="center", va="top", fontsize=10, color="#334155",
            linespacing=1.5)

    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    out = OUTPUT_DIR / "methods_overview_schematic.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
