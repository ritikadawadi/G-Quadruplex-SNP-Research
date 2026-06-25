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
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

REGIONS = {
    "5utr": {
        "label": "5'UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "12_top_genes_clinvar_qgrs_table.tsv",
        "header": "#2f7f5f",
        "accent": "#e6f4ea",
    },
    "3utr": {
        "label": "3'UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "12_top_genes_clinvar_qgrs_table.tsv",
        "header": "#6a51a3",
        "accent": "#efe8fb",
    },
}

TOP_N = 6


def first_terms(block: str, limit: int = 2, width: int = 26) -> str:
    if not isinstance(block, str) or not block.strip():
        return "—"
    lines = [line.strip("• ").strip() for line in block.splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith("(+") and "other terms" not in line]
    if not lines:
        return "—"
    lines = lines[:limit]
    wrapped = []
    for line in lines:
        if len(line) <= width:
            wrapped.append(line)
        else:
            wrapped.append(line[: width - 1] + "…")
    return "\n".join(wrapped)


def simplify_role(role: str, width: int = 16) -> str:
    role = str(role).strip()
    if not role or role.lower() == "nan":
        return "—"
    if len(role) <= width:
        return role
    return role[: width - 1] + "…"


def avg_delta(row: pd.Series) -> str:
    n = int(row["N_SNV"])
    if n == 0:
        return "0"
    mean = (float(row["SUM_POS"]) + float(row["SUM_NEG"])) / n
    return f"{mean:+.1f}"


def load_panel_df(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.head(TOP_N).copy()
    df["Key conditions"] = df["CLNDN"].map(first_terms)
    df["Role"] = df["COSMIC_ROLE"].map(simplify_role)
    df["Avg Δ"] = df.apply(avg_delta, axis=1)
    df["Σ|Δ|"] = df["TOTAL_ABS"].map(lambda x: f"{float(x):.1f}")
    df["SNVs"] = df["N_SNV"].astype(int).astype(str)
    return df[["GENE", "Role", "Key conditions", "Avg Δ", "Σ|Δ|", "SNVs"]]


def draw_panel(ax, df: pd.DataFrame, title: str, header_color: str, accent_color: str) -> None:
    ax.axis("off")
    ax.text(
        0.5,
        1.04,
        title,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
        color="#1f1f1f",
    )

    headers = ["Gene", "Role", "Key conditions", "Avg Δ G4\nper SNV", "Σ|Δ|", "SNVs"]
    cells = df.values.tolist()
    col_widths = [0.12, 0.15, 0.38, 0.12, 0.10, 0.10]
    tbl = ax.table(
        cellText=cells,
        colLabels=headers,
        cellLoc="left",
        colWidths=col_widths,
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 2.2)

    zebra_a = "#ffffff"
    zebra_b = "#f7f8fa"

    for (i, j), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d7dce2")
        cell.set_linewidth(0.7)
        if i == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(color="white", fontweight="bold", fontsize=9)
            cell.PAD = 0.03
            continue

        cell.set_facecolor(zebra_a if i % 2 else zebra_b)
        if j == 0:
            cell.set_text_props(fontweight="bold", fontsize=9.5)
        if j == 2:
            cell.set_facecolor(accent_color)
        if j in (4, 5):
            cell.set_text_props(ha="right", fontfamily="DejaVu Sans Mono", fontsize=9.5)
        if j == 3:
            txt = cell.get_text().get_text().strip()
            try:
                val = float(txt)
                if val < 0:
                    cell.set_text_props(color="#c05a00", fontweight="bold")
                elif val > 0:
                    cell.set_text_props(color="#127a3d", fontweight="bold")
            except ValueError:
                pass


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.8))
    fig.patch.set_facecolor("white")

    for ax, meta in zip(axes, REGIONS.values()):
        df = load_panel_df(meta["table"])
        draw_panel(ax, df, meta["label"], meta["header"], meta["accent"])

    fig.suptitle(
        "Top genes with the strongest G4-score disruption",
        fontsize=16,
        fontweight="bold",
        y=0.97,
    )
    fig.text(
        0.5,
        0.92,
        "Avg Δ G4 per SNV = mean change in max G-score per variant (negative = G4 weakened).  Σ|Δ| = total absolute change.",
        ha="center",
        fontsize=9,
        color="#555555",
    )
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.90])

    outfile = OUTPUT_DIR / "combined_top_genes_g4_disruption_5utr_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
