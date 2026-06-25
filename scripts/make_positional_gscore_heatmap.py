#!/usr/bin/env python3
"""
Gene-schematic figure with heatmap strips, dot clusters, chromosome
distribution, and region-length context showing where G-score-disrupting
SNVs concentrate along each genomic region.
"""
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

REGIONS = {
    "tss": {
        "table": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
        "len_col": "TSS_UPSTREAM500_LENGTH",
        "color": "#e6550d",
        "label": "Upstream",
        "short": "Upstream\n(500 bp)",
    },
    "5utr": {
        "table": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
        "len_col": "FIVE_UTR_LENGTH",
        "color": "#2ca25f",
        "label": "5\u2032UTR",
        "short": "5\u2032UTR",
    },
    "3utr": {
        "table": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
        "len_col": "THREE_UTR_LENGTH",
        "color": "#756bb1",
        "label": "3\u2032UTR",
        "short": "3\u2032UTR",
    },
}

HEATMAP_BINS = 40
DOT_CMAP = "YlOrRd"
HEATMAP_CMAP = "YlOrRd"
CDS_COLOR = "#aaaaaa"
SEGMENT_GAP = 0.015
# Gene schematic bar geometry (shared by draw_gene_schematic and Fig 1 axis labels)
GENE_SCHEMATIC_BAR_Y = 0.0
GENE_SCHEMATIC_BAR_H = 0.38
# Fig 1 auxiliary grey text (larger for poster readability)
FIG1_GREY_NOTE = 14
FIG1_GREY_MED = 12
FIG1_GREY_LANDMARK = 13
FIG1_GREY_CBAR_TICK = 11
FIG1_GENE_BLOCK_LABEL = 12
# Raise 0%/100% + double-arrow above UTR gene blocks (data coords)
FIG1_UTR_PCT_LIFT = 0.048
# UTR n= counts on dot panel (data coords, 0–1 y in draw_dot_cluster)
FIG1_DOTS_UTR_N_Y = 0.035
FIG1_UTR_0100_FS = 13
FIG1_HEATMAP_N_FS = 13


def compute_proportional_widths(transcripts: pd.DataFrame) -> dict:
    """Equal widths for the three main regions, with a smaller CDS spacer."""
    return {"tss": 1.0, "5utr": 1.0, "cds": 0.35, "3utr": 1.0}


def region_length_table(transcripts: pd.DataFrame) -> dict:
    """Return dict of {region: {median, iqr_lo, iqr_hi, pct_of_total}}."""
    fu = transcripts.loc[transcripts["FIVE_UTR_LENGTH"] > 0, "FIVE_UTR_LENGTH"]
    tu = transcripts.loc[transcripts["THREE_UTR_LENGTH"] > 0, "THREE_UTR_LENGTH"]
    med_up, med_5, med_3 = 500, int(fu.median()), int(tu.median())
    total = med_up + med_5 + med_3
    return {
        "tss": {"median": med_up, "iqr": "fixed", "pct": 100 * med_up / total,
                "range": f"500 nt (fixed)"},
        "5utr": {"median": med_5,
                  "iqr": f"{int(fu.quantile(0.25))}\u2013{int(fu.quantile(0.75))} nt",
                  "pct": 100 * med_5 / total,
                  "range": f"{int(fu.min())}\u2013{int(fu.max())} nt, med {med_5}"},
        "3utr": {"median": med_3,
                  "iqr": f"{int(tu.quantile(0.25))}\u2013{int(tu.quantile(0.75))} nt",
                  "pct": 100 * med_3 / total,
                  "range": f"{int(tu.min())}\u2013{int(tu.max())} nt, med {med_3}"},
    }


def load_and_normalise(region_key: str, meta: dict, transcripts: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(meta["table"], sep="\t")
    df = df.merge(
        transcripts[["TRANSCRIPT_ID", meta["len_col"]]],
        left_on="TRANSCRIPT",
        right_on="TRANSCRIPT_ID",
        how="left",
    )
    region_len = df[meta["len_col"]]
    df["frac_pos"] = df["TRANSCRIPT_POS1"] / region_len
    df["frac_pos"] = df["frac_pos"].clip(0, 1)
    df["abs_delta"] = df["DELTA_MAX_GSCORE"].abs()
    df["changed"] = df["abs_delta"] > 0
    return df


def build_seg_bounds(seg_widths: dict, gap: float = SEGMENT_GAP) -> dict:
    seg_order = ["tss", "5utr", "cds", "3utr"]
    total_gap = gap * (len(seg_order) - 1)
    total_w = sum(seg_widths[s] for s in seg_order) + total_gap
    scale = 1.0 / total_w
    bounds = {}
    cursor = 0.0
    for s in seg_order:
        w = seg_widths[s] * scale
        bounds[s] = (cursor, cursor + w)
        cursor += w + gap * scale
    return bounds


def draw_gene_schematic(ax, seg_bounds: dict, len_tbl: dict, show_annotations: bool = True):
    """Gene schematic with proportional blocks and optional median length / %total annotations."""
    bar_y, bar_h = GENE_SCHEMATIC_BAR_Y, GENE_SCHEMATIC_BAR_H
    for key in ["tss", "5utr", "cds", "3utr"]:
        x0, x1 = seg_bounds[key]
        if key == "cds":
            color, label = CDS_COLOR, "CDS"
        else:
            color = REGIONS[key]["color"]
            label = REGIONS[key]["short"]
        rect = mpatches.FancyBboxPatch(
            (x0, bar_y), x1 - x0, bar_h,
            boxstyle="round,pad=0.005",
            facecolor=color, edgecolor="white", linewidth=1.2, zorder=3,
        )
        ax.add_patch(rect)
        ax.text((x0 + x1) / 2, bar_y + bar_h / 2, label,
                ha="center", va="center", fontsize=FIG1_GENE_BLOCK_LABEL, fontweight="bold",
                color="white", zorder=4)

    if show_annotations:
        for key in ["tss", "5utr", "3utr"]:
            x0, x1 = seg_bounds[key]
            info = len_tbl[key]
            pct_txt = f"{info['pct']:.0f}% of non-coding span"
            ax.text((x0 + x1) / 2, bar_y + bar_h + 0.04, info["range"],
                    ha="center", va="bottom", fontsize=FIG1_GREY_LANDMARK, color="#444444")
            ax.text((x0 + x1) / 2, bar_y + bar_h + 0.16, pct_txt,
                    ha="center", va="bottom", fontsize=FIG1_GREY_LANDMARK, color="#666666", fontstyle="italic")

    cx0, cx1 = seg_bounds["cds"]
    cds_w = cx1 - cx0
    tss_x0, tss_x1 = seg_bounds["tss"]
    ax.text(
        (tss_x0 + tss_x1) / 2, bar_y - 0.07,
        "\u2190 500 bp \u2192",
        ha="center", va="top", fontsize=FIG1_GREY_LANDMARK, color="#555555",
    )
    landmarks = [
        (seg_bounds["tss"][1], "TSS", "center"),
        (cx0 + 0.18 * cds_w, "Start\ncodon", "center"),
        (cx1 - 0.18 * cds_w, "Stop\ncodon", "center"),
        (seg_bounds["3utr"][1], "Poly(A) \u2192", "right"),
    ]
    for lx, ltxt, ha in landmarks:
        ax.text(lx, bar_y - 0.07, ltxt, ha=ha, va="top", fontsize=FIG1_GREY_LANDMARK, color="#555555")

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.35, bar_h + 0.35)
    ax.axis("off")


def draw_heatmap_strip(ax, frac_positions, seg_x0, seg_x1, n_bins=HEATMAP_BINS):
    counts, _ = np.histogram(frac_positions, bins=n_bins, range=(0, 1))
    seg_w = seg_x1 - seg_x0
    bin_w = seg_w / n_bins
    cmap = plt.get_cmap(HEATMAP_CMAP)
    norm = Normalize(vmin=0, vmax=max(counts.max(), 1))
    for i, c in enumerate(counts):
        rect = mpatches.Rectangle(
            (seg_x0 + i * bin_w, 0.0), bin_w, 1.0,
            facecolor=cmap(norm(c)), edgecolor="none",
        )
        ax.add_patch(rect)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.1, 1.1)
    ax.axis("off")
    return norm


def genomic_variant_span_mb(rd: pd.DataFrame) -> str:
    """Min–max variant genomic position in Mb (same summary as region stats table)."""
    if rd is None or len(rd) == 0:
        return "\u2014"
    lo = float(rd["VARIANT_POS_1BASED"].min()) / 1e6
    hi = float(rd["VARIANT_POS_1BASED"].max()) / 1e6
    return f"{lo:.1f}\u2013{hi:.1f} Mb"


def save_positional_distribution_figure(
    region_data: dict,
    seg_bounds: dict,
    len_tbl: dict,
    vmin: float,
    vmax: float,
    outfile: Path,
    *,
    add_genomic_span_labels: bool = False,
) -> None:
    """Dots + SNV density + gene schematic; optional genomic span (Mb) above each region."""
    fig1 = plt.figure(figsize=(16, 7), facecolor="white")
    gs1 = gridspec.GridSpec(
        3, 1, figure=fig1, height_ratios=[3, 0.5, 1.2], hspace=0.06,
        left=0.05, right=0.88, top=0.85, bottom=0.10,
    )
    ax_dots = fig1.add_subplot(gs1[0])
    ax_heat = fig1.add_subplot(gs1[1])
    ax_gene = fig1.add_subplot(gs1[2])

    sc = None
    for key in ["tss", "5utr", "3utr"]:
        rd = region_data[key]
        ret = draw_dot_cluster(
            ax_dots, rd["frac_pos"].values, rd["abs_delta"].values,
            seg_bounds[key][0], seg_bounds[key][1], vmin, vmax,
        )
        if ret is not None:
            sc = ret
    for key in ["tss", "5utr", "3utr"]:
        x0, x1 = seg_bounds[key]
        ax_dots.axvline(x0, color="#cccccc", lw=0.5, ls="--", zorder=1)
        ax_dots.axvline(x1, color="#cccccc", lw=0.5, ls="--", zorder=1)
    n_total = sum(len(rd) for rd in region_data.values())
    ax_dots.text(
        0.99, 0.97,
        f"n = {n_total} SNVs with \u0394 max G-score \u2260 0",
        transform=ax_dots.transAxes, ha="right", va="top",
        fontsize=FIG1_GREY_NOTE, color="#333333", fontweight="semibold",
    )
    for key in ("tss", "5utr", "3utr"):
        x0, x1 = seg_bounds[key]
        n_reg = len(region_data[key])
        ax_dots.text(
            (x0 + x1) / 2, FIG1_DOTS_UTR_N_Y, f"n = {n_reg}",
            ha="center", va="bottom", fontsize=FIG1_HEATMAP_N_FS,
            color="#2d2d2d", fontweight="semibold", zorder=6,
        )
    # n= labels moved to dot panel; nothing printed under heatmap

    for key in ["tss", "5utr", "3utr"]:
        rd = region_data[key]
        draw_heatmap_strip(ax_heat, rd["frac_pos"].values,
                           seg_bounds[key][0], seg_bounds[key][1])
    cx0, cx1 = seg_bounds["cds"]
    ax_heat.add_patch(mpatches.Rectangle(
        (cx0, 0.0), cx1 - cx0, 1.0, facecolor="#eeeeee", edgecolor="none",
    ))
    ax_heat.text(0.0, 0.5, "SNV\ndensity", ha="right", va="center", fontsize=FIG1_GREY_MED,
                 color="#555555", fontweight="semibold", transform=ax_heat.transAxes)
    ax_heat.set_ylim(-0.1, 1.06)

    draw_gene_schematic(ax_gene, seg_bounds, len_tbl, show_annotations=False)
    _by, _bh = GENE_SCHEMATIC_BAR_Y, GENE_SCHEMATIC_BAR_H
    if add_genomic_span_labels:
        _span_y = _by + _bh + 0.055
        for key in ["tss", "5utr", "3utr"]:
            x0, x1 = seg_bounds[key]
            ax_gene.text(
                (x0 + x1) / 2, _span_y,
                genomic_variant_span_mb(region_data[key]),
                ha="center", va="bottom", fontsize=FIG1_GREY_MED,
                color="#2d2d2d", fontweight="semibold", zorder=5,
            )
        _frac_y = _by + _bh + 0.13
        ax_gene.set_ylim(-0.35, _bh + 0.52)
    else:
        _frac_y = _by + _bh + 0.03

    _utr_pct_y = _frac_y + FIG1_UTR_PCT_LIFT
    _utr_arrow_pad = 0.042
    _len_label_dy = 0.042          # gap: arrow line → "length" label above it
    _length_word = "length"

    # Upstream arrow and length label
    tss_x0, tss_x1 = seg_bounds["tss"]
    _tss_arrow_pad = 0.015
    ax_gene.add_patch(mpatches.FancyArrowPatch(
        (tss_x0 + _tss_arrow_pad, _utr_pct_y),
        (tss_x1 - _tss_arrow_pad, _utr_pct_y),
        arrowstyle="<->", mutation_scale=14,
        color="#555555", linewidth=1.45, zorder=4,
    ))
    ax_gene.text(
        (tss_x0 + tss_x1) / 2, _utr_pct_y + _len_label_dy,
        _length_word,
        ha="center", va="bottom", fontsize=FIG1_UTR_0100_FS,
        color="#333333", fontweight="bold", zorder=5,
    )

    # 5'UTR and 3'UTR arrows, 0%/100% labels, and length labels
    for key in ("5utr", "3utr"):
        x0, x1 = seg_bounds[key]
        ax_gene.add_patch(mpatches.FancyArrowPatch(
            (x0 + _utr_arrow_pad, _utr_pct_y),
            (x1 - _utr_arrow_pad, _utr_pct_y),
            arrowstyle="<->", mutation_scale=14,
            color="#555555", linewidth=1.45, zorder=4,
        ))
        ax_gene.text(
            x0 + 0.006, _utr_pct_y, "0%",
            fontsize=FIG1_UTR_0100_FS, color="#333333", va="center", ha="left",
            zorder=5, fontweight="bold",
        )
        ax_gene.text(
            x1 - 0.006, _utr_pct_y, "100%",
            fontsize=FIG1_UTR_0100_FS, color="#333333", va="center", ha="right",
            zorder=5, fontweight="bold",
        )
        ax_gene.text(
            (x0 + x1) / 2, _utr_pct_y + _len_label_dy,
            _length_word,
            ha="center", va="bottom", fontsize=FIG1_UTR_0100_FS,
            color="#333333", fontweight="bold", zorder=5,
        )
    _gene_ylim_top = max(float(ax_gene.get_ylim()[1]), _utr_pct_y + _len_label_dy + 0.08)
    ax_gene.set_ylim(ax_gene.get_ylim()[0], _gene_ylim_top)

    if sc is not None:
        cbar_ax = fig1.add_axes([0.91, 0.38, 0.012, 0.40])
        sm = ScalarMappable(cmap=DOT_CMAP, norm=Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = fig1.colorbar(sm, cax=cbar_ax)
        cbar.set_label("|\u0394 max G-score|", fontsize=FIG1_GREY_MED, color="#222222")
        cbar.ax.tick_params(labelsize=FIG1_GREY_CBAR_TICK, colors="#333333")

    fig1.text(
        0.46, 0.98,
        "Positional distribution of G-score-disrupting SNVs across gene regions",
        ha="center", va="top", fontsize=24, fontweight="bold",
    )
    fig1.text(
        0.46, 0.898,
        "Where in each region do G-score-disrupting SNVs fall?",
        ha="center", va="top", fontsize=14, fontweight="bold", color="#333333",
    )

    footer = (
        "Each dot is one SNV whose strongest predicted G-score changed (|\u0394| > 0). "
        "Heatmap strip shows binned density. "
        "Position is normalised within each region (0 = region start, 1 = region end)."
    )
    if add_genomic_span_labels:
        footer += (
            " Values above Upstream and UTR blocks: min\u2013max genomic coordinate span (Mb) "
            "of those SNVs (positions are not necessarily contiguous on one chromosome)."
        )
    fig1.text(
        0.46, 0.02,
        footer,
        ha="center", va="bottom", fontsize=FIG1_GREY_NOTE, color="#3d3d3d",
    )

    fig1.savefig(outfile, dpi=300, bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig1)
    print(f"Wrote {outfile}")


def draw_dot_cluster(ax, frac_positions, abs_deltas, seg_x0, seg_x1, vmin, vmax):
    if len(frac_positions) == 0:
        return None
    seg_w = seg_x1 - seg_x0
    mapped_x = seg_x0 + frac_positions * seg_w
    rng = np.random.default_rng(42)
    jitter_y = rng.uniform(0.08, 0.92, size=len(mapped_x))
    sc = ax.scatter(
        mapped_x, jitter_y,
        c=abs_deltas, cmap=DOT_CMAP, norm=Normalize(vmin=vmin, vmax=vmax),
        s=18, alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3,
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.12, 1.05)
    ax.axis("off")
    return sc


def draw_chrom_bar(ax, region_data: dict):
    """Small grouped bar chart: chromosome distribution of disrupted SNVs across regions."""
    all_chroms = set()
    for rd in region_data.values():
        all_chroms.update(rd["VARIANT_CHROM"].astype(str).unique())

    def chrom_sort_key(c):
        try:
            return (0, int(c))
        except ValueError:
            return (1, c)
    chroms = sorted(all_chroms, key=chrom_sort_key)

    x = np.arange(len(chroms))
    width = 0.28
    region_order = ["tss", "5utr", "3utr"]
    for idx, key in enumerate(region_order):
        counts = []
        for c in chroms:
            counts.append((region_data[key]["VARIANT_CHROM"].astype(str) == c).sum())
        offset = (idx - 1) * width
        ax.bar(x + offset, counts, width=width, color=REGIONS[key]["color"],
               label=REGIONS[key]["label"], edgecolor="white", linewidth=0.4)
    ax.set_xticks(x, [f"chr{c}" for c in chroms], fontsize=6, rotation=60, ha="right")
    ax.set_ylabel("Disrupted SNVs", fontsize=8)
    ax.set_title("Chromosomal distribution of G-score-disrupting SNVs", fontsize=18, fontweight="bold", pad=8)
    ax.legend(fontsize=7, frameon=True, loc="upper right", ncol=1)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlim(-0.6, len(chroms) - 0.4)
    ax.grid(axis="y", alpha=0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_stats_box(ax, region_data: dict, len_tbl: dict):
    """Summary statistics table in a dedicated axes."""
    ax.axis("off")
    header = ["Region", "Disrupted\nSNVs", "Unique\ngenes", "Median\nlength",
              "% of non-\ncoding span", "Mean |Δ|",
              "Genomic pos\nrange (Mb)"]
    rows = []
    for key in ["tss", "5utr", "3utr"]:
        rd = region_data[key]
        n_genes = rd["GENE"].nunique()
        mean_delta = rd["abs_delta"].mean()
        info = len_tbl[key]
        pos_min = rd["VARIANT_POS_1BASED"].min() / 1e6
        pos_max = rd["VARIANT_POS_1BASED"].max() / 1e6
        rows.append([
            REGIONS[key]["label"],
            str(len(rd)),
            str(n_genes),
            f"{info['median']} nt",
            f"{info['pct']:.1f}%",
            f"{mean_delta:.1f}",
            f"{pos_min:.1f}\u2013{pos_max:.1f}",
        ])
    ax.set_title("Region summary", fontsize=18, fontweight="bold", pad=12)
    tbl = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.0, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#f0f0f0")
            cell.set_text_props(fontweight="bold")
        elif c == 0:
            cell.set_facecolor(REGIONS[["tss", "5utr", "3utr"][r - 1]]["color"] + "22")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    transcripts = pd.read_csv(
        PROJECT_ROOT / "data" / "regions" / "selected_transcripts.tsv", sep="\t",
    )

    region_all = {}
    region_data = {}
    for key, meta in REGIONS.items():
        df = load_and_normalise(key, meta, transcripts)
        region_all[key] = df
        region_data[key] = df[df["changed"]].copy()

    seg_widths = compute_proportional_widths(transcripts)
    seg_bounds = build_seg_bounds(seg_widths)
    len_tbl = region_length_table(transcripts)

    all_deltas = np.concatenate([rd["abs_delta"].values for rd in region_data.values()])
    vmin, vmax = 0, max(all_deltas.max(), 1) if len(all_deltas) else 1

    # ================================================================
    # FIGURE 1: Positional distribution (dots + heatmap + gene schematic)
    # ================================================================
    save_positional_distribution_figure(
        region_data, seg_bounds, len_tbl, vmin, vmax,
        OUTPUT_DIR / "positional_gscore_disruption_gene_schematic.png",
        add_genomic_span_labels=False,
    )
    save_positional_distribution_figure(
        region_data, seg_bounds, len_tbl, vmin, vmax,
        OUTPUT_DIR / "positional_gscore_disruption_gene_schematic_genomic_span.png",
        add_genomic_span_labels=True,
    )

    # ================================================================
    # FIGURE 2: Chromosomal distribution + region summary table (stacked)
    # ================================================================
    fig2 = plt.figure(figsize=(14, 9), facecolor="white")
    gs2 = gridspec.GridSpec(
        2, 1, figure=fig2, height_ratios=[2.0, 1.05], hspace=0.32,
        left=0.08, right=0.96, top=0.94, bottom=0.08,
    )
    ax_chrom = fig2.add_subplot(gs2[0])
    ax_stats = fig2.add_subplot(gs2[1])

    draw_chrom_bar(ax_chrom, region_data)
    draw_stats_box(ax_stats, region_data, len_tbl)

    outfile2 = OUTPUT_DIR / "chromosomal_distribution_region_summary.png"
    fig2.savefig(outfile2, dpi=300, bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig2)
    print(f"Wrote {outfile2}")


if __name__ == "__main__":
    main()
