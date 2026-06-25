#!/usr/bin/env python3
"""
Combined annotated bar chart: net G4-score change per gene with
COSMIC role + key disease annotation, for 5'UTR and 3'UTR.
"""
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

SUMMARY = {
    "5utr": PROJECT_ROOT / "results_compare_figures/5utr/tables/5utr_comparison_summary_copy.tsv",
    "3utr": PROJECT_ROOT / "results_compare_figures/3utr/tables/3utr_comparison_summary_copy.tsv",
}
CLINVAR_TABLE = {
    "5utr": PROJECT_ROOT / "results_compare_figures/5utr/tables/12_top_genes_clinvar_qgrs_table.tsv",
    "3utr": PROJECT_ROOT / "results_compare_figures/3utr/tables/12_top_genes_clinvar_qgrs_table.tsv",
}

N_LOSS = 4
N_GAIN = 2
CLR_LOSS = "#d7301f"
CLR_GAIN = "#2171b5"


def first_disease(clndn_block: str) -> str:
    if not isinstance(clndn_block, str) or not clndn_block.strip():
        return ""
    lines = [l.strip("• ").strip() for l in clndn_block.splitlines() if l.strip()]
    lines = [l for l in lines if not l.startswith("(+") and "other terms" not in l]
    if not lines:
        return ""
    term = lines[0]
    if len(term) > 28:
        term = term[:27] + "…"
    return term


def build_data(region: str) -> pd.DataFrame:
    summary = pd.read_csv(SUMMARY[region], sep="\t")
    g = summary.groupby("GENE").agg(
        net_delta=("DELTA_MAX_GSCORE", "sum"),
        n=("DELTA_MAX_GSCORE", "count"),
    ).reset_index()
    g["abs_net"] = g["net_delta"].abs()

    losses = g[g["net_delta"] < 0].sort_values("abs_net", ascending=False).head(N_LOSS)
    gains = g[g["net_delta"] > 0].sort_values("net_delta", ascending=False).head(N_GAIN)
    combo = pd.concat([losses, gains]).sort_values("net_delta")

    info = pd.read_csv(CLINVAR_TABLE[region], sep="\t")
    role_map = dict(zip(info["GENE"], info["COSMIC_ROLE"].fillna("")))
    disease_map = dict(zip(info["GENE"], info["CLNDN"].fillna("").map(first_disease)))

    combo["role"] = combo["GENE"].map(lambda g: role_map.get(g, ""))
    combo["disease"] = combo["GENE"].map(lambda g: disease_map.get(g, ""))
    return combo.reset_index(drop=True)


def draw_panel(ax, df: pd.DataFrame, title: str):
    y = np.arange(len(df))
    vals = df["net_delta"].values.astype(float)
    colors = [CLR_LOSS if v < 0 else CLR_GAIN for v in vals]

    ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.8, height=0.6)
    ax.axvline(0, color="#333333", linewidth=0.9)

    mag = max(abs(vals.min()), abs(vals.max())) if len(vals) else 1.0
    xlim = mag * 1.45
    ax.set_xlim(-xlim, xlim)

    gene_labels = []
    for _, row in df.iterrows():
        role = row["role"]
        label = row["GENE"]
        if role:
            label += f"  ({role})"
        gene_labels.append(label)

    ax.set_yticks(y, gene_labels)
    ax.tick_params(axis="y", labelsize=9)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=12, fontweight="bold", pad=6)
    ax.set_xlabel("Net Δ max G-score", fontsize=9)
    ax.grid(axis="x", alpha=0.2)

    for i, row in df.iterrows():
        val = float(row["net_delta"])
        n_snvs = int(row["n"])
        disease = row["disease"]

        offset = xlim * 0.02
        if val < 0:
            n_x, n_ha = val - offset, "right"
            d_x, d_ha = offset * 2, "left"
        else:
            n_x, n_ha = val + offset, "left"
            d_x, d_ha = -offset * 2, "right"

        ax.text(n_x, y[i], f"n={n_snvs}", ha=n_ha, va="center",
                fontsize=7.5, color="#555555")

        if disease:
            ax.text(d_x, y[i] + 0.28, disease, ha=d_ha, va="top",
                    fontsize=6.5, color="#777777", style="italic")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    data_5 = build_data("5utr")
    data_3 = build_data("3utr")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    fig.patch.set_facecolor("white")

    draw_panel(ax1, data_5, "5\u2032UTR")
    draw_panel(ax2, data_3, "3\u2032UTR")

    fig.suptitle(
        "Gene-level G4-score impact: losses vs gains",
        fontsize=14, fontweight="bold", y=1.01,
    )

    from matplotlib.patches import Patch
    legend = [
        Patch(facecolor=CLR_LOSS, label="G4 loss (score decrease)"),
        Patch(facecolor=CLR_GAIN, label="G4 gain (score increase)"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=2,
               fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.text(0.5, -0.10,
             "Bars = net Δ max G-score summed across SNVs.  "
             "Gene labels include COSMIC role.  n = variant count.",
             ha="center", fontsize=8, color="#777777")

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    out = OUTPUT_DIR / "gene_g4_overview_5utr_3utr.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
