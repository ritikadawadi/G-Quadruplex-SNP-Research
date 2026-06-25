#!/usr/bin/env python3
"""
Poster tables for top genes (by sum |Δ max g-score|), with separate COSMIC vs
ClinVar columns and readable styling. See also notebooks/top_gene_clinvar_tables.ipynb.

Outputs:
  - results_compare_figures/{region}/figures/12_top_genes_clinvar_qgrs_table.png
  - results_compare_figures/{region}/tables/12_top_genes_clinvar_qgrs_table.tsv
"""
from __future__ import annotations

import os
import re
from collections import Counter
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)

import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONTEXT_ROOT = PROJECT_ROOT / "results_cancer_context_figures"
COMPARE_FIG_ROOT = PROJECT_ROOT / "results_compare_figures"
FASTA_ROOT = PROJECT_ROOT / "initial_fasta"

REGIONS = {
    "tss": {
        "label": "TSS (upstream 500 bp)",
        "annotated": CONTEXT_ROOT / "tss" / "tables" / "tss_cosmic_annotated_comparison.tsv",
        "qgrs_hits": FASTA_ROOT / "results_tss" / "tss_upstream500_reference.qgrs_hits.tsv",
    },
    "5utr": {
        "label": "5′UTR",
        "annotated": CONTEXT_ROOT / "5utr" / "tables" / "5utr_cosmic_annotated_comparison.tsv",
        "qgrs_hits": FASTA_ROOT / "results_5utr" / "5utr_reference.qgrs_hits.tsv",
    },
    "3utr": {
        "label": "3′UTR",
        "annotated": CONTEXT_ROOT / "3utr" / "tables" / "3utr_cosmic_annotated_comparison.tsv",
        "qgrs_hits": FASTA_ROOT / "results_3utr" / "3utr_reference.qgrs_hits.tsv",
    },
}

TOP_N = 15


def split_clndn(raw: str) -> list[str]:
    if not isinstance(raw, str) or not raw.strip():
        return []
    parts = re.split(r"[|;,]", raw)
    out = []
    for p in parts:
        t = p.replace("_", " ").strip()
        if t and t.lower() not in {"not provided", "not specified", "nan"}:
            out.append(t)
    return out


def parse_overlap_ids(val) -> list[str]:
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s:
        return []
    return [p for p in re.split(r"[,;\s]+", s) if p.strip().isdigit()]


def clinsig_compact(series: pd.Series, max_cats: int = 3) -> str:
    vals = series.dropna().astype(str).str.strip().replace("", pd.NA).dropna()
    vals = vals[vals.str.lower() != "nan"]
    if vals.empty:
        return "—"
    cnt = vals.value_counts()
    parts = []
    for sig, c in cnt.head(max_cats).items():
        short = sig.replace("_", " ")[:40] + ("…" if len(sig) > 40 else "")
        parts.append(f"{short} ({c})")
    if len(cnt) > max_cats:
        parts.append(f"+{len(cnt) - max_cats} more")
    return "\n".join(parts)


def diseases_compact(series: pd.Series, top_n: int = 5) -> str:
    bag: Counter[str] = Counter()
    for raw in series.dropna():
        for d in split_clndn(str(raw)):
            bag[d] += 1
    if not bag:
        return "—"
    lines = []
    for name, _ in bag.most_common(top_n):
        t = name if len(name) <= 52 else name[:49] + "…"
        lines.append(f"• {t}")
    if len(bag) > top_n:
        lines.append(f"• (+{len(bag) - top_n} other terms)")
    return "\n".join(lines)


def build_top_gene_df(region_key: str) -> pd.DataFrame:
    meta = REGIONS[region_key]
    df = pd.read_csv(meta["annotated"], sep="\t")
    df["DELTA_MAX_GSCORE"] = pd.to_numeric(df["DELTA_MAX_GSCORE"], errors="coerce").fillna(0)

    gene_rank = (
        df.groupby("GENE")
        .agg(total_abs=("DELTA_MAX_GSCORE", lambda s: s.abs().sum()))
        .sort_values("total_abs", ascending=False)
        .head(TOP_N)
    )
    top_genes = gene_rank.index.tolist()

    qh = pd.read_csv(meta["qgrs_hits"], sep="\t")
    qgrs_per_gene = qh.groupby("GENE").size().reindex(top_genes).fillna(0).astype(int)

    rows = []
    for gene in top_genes:
        g = df[df["GENE"] == gene]
        d = g["DELTA_MAX_GSCORE"]
        sum_pos = float(d[d > 0].sum())
        sum_neg = float(d[d < 0].sum())
        overlap_ids: set[str] = set()
        for v in g["QGRS_OVERLAP_IDS"]:
            overlap_ids.update(parse_overlap_ids(v))
        first = g.iloc[0]
        cosmic = (
            str(first["ROLE_IN_CANCER"]).strip()
            if "ROLE_IN_CANCER" in g.columns and pd.notna(first.get("ROLE_IN_CANCER"))
            else "—"
        )
        rows.append(
            {
                "GENE": gene,
                "COSMIC_ROLE": cosmic,
                "CLINSIG": clinsig_compact(g["CLNSIG"]) if "CLNSIG" in g.columns else "—",
                "CLNDN": diseases_compact(g["CLNDN"]) if "CLNDN" in g.columns else "—",
                "SUM_POS": sum_pos,
                "SUM_NEG": sum_neg,
                "N_SNV": len(g),
                "REF_QGRS": int(qgrs_per_gene.get(gene, 0)),
                "QGRS_WITH_SNP": len(overlap_ids),
                "TOTAL_ABS": float(gene_rank.loc[gene, "total_abs"]),
            }
        )
    return pd.DataFrame(rows)


def draw_clean_table(region_key: str, region_label: str, out_df: pd.DataFrame) -> None:
    fig_dir = COMPARE_FIG_ROOT / region_key / "figures"
    tbl_dir = COMPARE_FIG_ROOT / region_key / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    headers = [
        "Gene",
        "Cancer role\n(COSMIC)",
        "Clinical significance\n(ClinVar CLNSIG)",
        "Associated conditions\n(ClinVar CLNDN)",
        "Σ Δ max\ng-score\n(+)",
        "Σ Δ max\ng-score\n(−)",
        "Ref.\nQGRS",
        "QGRS with\nG-SNP",
    ]

    cells = []
    for _, r in out_df.iterrows():
        cells.append(
            [
                r["GENE"],
                r["COSMIC_ROLE"],
                r["CLINSIG"],
                r["CLNDN"],
                f"{r['SUM_POS']:.1f}" if r["SUM_POS"] else "0",
                f"{r['SUM_NEG']:.1f}" if r["SUM_NEG"] else "0",
                str(r["REF_QGRS"]),
                str(r["QGRS_WITH_SNP"]),
            ]
        )

    fig = plt.figure(figsize=(26, 12), constrained_layout=False)
    fig.patch.set_facecolor("#f8f9fa")
    ax = fig.add_axes([0.04, 0.14, 0.92, 0.78])
    ax.axis("off")
    fig.suptitle(
        f"{region_label}: top {TOP_N} genes by Σ|Δ max g-score|",
        fontsize=18,
        fontweight="bold",
        y=0.97,
        color="#1a1a1a",
    )
    ax.text(
        0.5,
        1.02,
        "G-related SNVs in COSMIC cancer genes — directional change in QGRS max g-score vs wild-type",
        transform=ax.transAxes,
        ha="center",
        fontsize=11,
        color="#444",
    )

    col_widths = [0.07, 0.11, 0.15, 0.35, 0.08, 0.08, 0.08, 0.08]
    tbl = ax.table(
        cellText=cells,
        colLabels=headers,
        loc="center",
        cellLoc="left",
        colWidths=col_widths,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 2.85)

    hdr_color = "#2c3e50"
    zebra_a = "#ffffff"
    zebra_b = "#f0f4f8"
    green = "#b8e0b8"
    orange = "#f4c090"

    for (i, j), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cfd8dc")
        cell.set_linewidth(0.8)
        if i == 0:
            cell.set_facecolor(hdr_color)
            cell.set_text_props(color="white", fontweight="bold", fontsize=10)
            if j in (4, 5, 6, 7):
                cell.get_text().set_ha("right")
            cell.set_height(0.055)
            cell.PAD = 0.04
            continue
        row_bg = zebra_a if i % 2 else zebra_b
        cell.set_facecolor(row_bg)
        if j == 0:
            cell.set_text_props(fontweight="bold", fontsize=10)
        if j in (4, 5, 6, 7):
            cell.set_text_props(fontfamily="DejaVu Sans Mono", fontsize=10, ha="right")
        if j == 4:
            cell.set_facecolor(green)
        elif j == 5:
            cell.set_facecolor(orange)

    foot = (
        "Σ Δ max g-score = sum of (mutant − wild-type) max g-score over G-related SNVs for that gene in this region.\n"
        "(+) column = sum of positive deltas only; (−) column = sum of negative deltas. "
        "Ref. QGRS = wild-type motif count in region FASTA. "
        "QGRS with G-SNP = distinct QGRS IDs overlapping ≥1 variant in the comparison table."
    )
    fig.text(0.5, 0.04, foot, ha="center", va="center", fontsize=10, color="#333", linespacing=1.35)

    fig.savefig(
        fig_dir / "12_top_genes_clinvar_qgrs_table.png",
        dpi=220,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )
    plt.close(fig)

    tsv_cols = [
        "GENE",
        "COSMIC_ROLE",
        "CLINSIG",
        "CLNDN",
        "SUM_POS",
        "SUM_NEG",
        "N_SNV",
        "REF_QGRS",
        "QGRS_WITH_SNP",
        "TOTAL_ABS",
    ]
    out_df[tsv_cols].to_csv(tbl_dir / "12_top_genes_clinvar_qgrs_table.tsv", sep="\t", index=False)


def main():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    for key, meta in REGIONS.items():
        if not meta["annotated"].is_file():
            raise FileNotFoundError(meta["annotated"])
        if not meta["qgrs_hits"].is_file():
            raise FileNotFoundError(meta["qgrs_hits"])
        out_df = build_top_gene_df(key)
        draw_clean_table(key, meta["label"], out_df)
        print(f"Wrote {COMPARE_FIG_ROOT / key / 'figures' / '12_top_genes_clinvar_qgrs_table.png'}")


if __name__ == "__main__":
    main()
