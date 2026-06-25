#!/usr/bin/env python3
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/Users/ritika/Desktop/pilot/results_compare_figures/.matplotlib")

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path("/Users/ritika/Desktop/pilot")
COMPARE_ROOT = PROJECT_ROOT / "compare_mutant_wild"
OUTPUT_ROOT = PROJECT_ROOT / "results_compare_figures"
MPL_DIR = OUTPUT_ROOT / ".matplotlib"

REGIONS = {
    "tss": {
        "label": "TSS",
        "summary": COMPARE_ROOT / "tss" / "tss_upstream500.comparison_summary.tsv",
        "overall": COMPARE_ROOT / "tss" / "tss_upstream500.overall_summary.tsv",
    },
    "5utr": {
        "label": "5'UTR",
        "summary": COMPARE_ROOT / "5utr" / "5utr.comparison_summary.tsv",
        "overall": COMPARE_ROOT / "5utr" / "5utr.overall_summary.tsv",
    },
    "3utr": {
        "label": "3'UTR",
        "summary": COMPARE_ROOT / "3utr" / "3utr.comparison_summary.tsv",
        "overall": COMPARE_ROOT / "3utr" / "3utr.overall_summary.tsv",
    },
}


def ensure_dirs():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    MPL_DIR.mkdir(parents=True, exist_ok=True)
    for key in REGIONS:
        region_dir = OUTPUT_ROOT / key
        (region_dir / "figures").mkdir(parents=True, exist_ok=True)
        (region_dir / "tables").mkdir(parents=True, exist_ok=True)


def configure_plotting():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        plt.style.use("ggplot")
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["axes.labelsize"] = 12
    plt.rcParams["xtick.labelsize"] = 10
    plt.rcParams["ytick.labelsize"] = 10


def load_region_df(region_key):
    df = pd.read_csv(REGIONS[region_key]["summary"], sep="\t")
    df["DELTA_QGRS_COUNT"] = pd.to_numeric(df["DELTA_QGRS_COUNT"])
    df["DELTA_MAX_GSCORE"] = pd.to_numeric(df["DELTA_MAX_GSCORE"])
    df["REFERENCE_QGRS_COUNT"] = pd.to_numeric(df["REFERENCE_QGRS_COUNT"])
    df["MUTANT_QGRS_COUNT"] = pd.to_numeric(df["MUTANT_QGRS_COUNT"])
    df["REFERENCE_MAX_GSCORE"] = pd.to_numeric(df["REFERENCE_MAX_GSCORE"])
    df["MUTANT_MAX_GSCORE"] = pd.to_numeric(df["MUTANT_MAX_GSCORE"])
    df["REFERENCE_MAX_TETRADS"] = pd.to_numeric(df["REFERENCE_MAX_TETRADS"])
    df["MUTANT_MAX_TETRADS"] = pd.to_numeric(df["MUTANT_MAX_TETRADS"])
    df["CLNSIG_PRIMARY"] = (
        df["CLNSIG"]
        .fillna("")
        .astype(str)
        .str.split("|").str[0]
        .str.split(",").str[0]
        .replace("", "not_reported")
    )
    return df


def save_table(df, path):
    df.to_csv(path, sep="\t", index=False)


def save_series(series, path, index_name, value_name):
    df = series.rename(value_name).reset_index().rename(columns={"index": index_name})
    df.to_csv(path, sep="\t", index=False)


def make_overview_figure(df, region_key, region_label, out_dir):
    effect_counts = df["EFFECT_CLASS"].value_counts()
    qgrs_change = pd.Series(
        {
            "More QGRS": (df["DELTA_QGRS_COUNT"] > 0).sum(),
            "Fewer QGRS": (df["DELTA_QGRS_COUNT"] < 0).sum(),
            "Same QGRS": (df["DELTA_QGRS_COUNT"] == 0).sum(),
        }
    )
    gscore_change = pd.Series(
        {
            "Higher Max G-score": (df["DELTA_MAX_GSCORE"] > 0).sum(),
            "Lower Max G-score": (df["DELTA_MAX_GSCORE"] < 0).sum(),
            "Same Max G-score": (df["DELTA_MAX_GSCORE"] == 0).sum(),
        }
    )
    motif_flags = pd.Series(
        {
            "Top Sequence Changed": (df["TOP_SEQUENCE_CHANGED"] == "y").sum(),
            "Top Sequence Unchanged": (df["TOP_SEQUENCE_CHANGED"] == "n").sum(),
            "Gained Motif Instances": (df["GAINED_MOTIF_INSTANCES"] > 0).sum(),
            "Lost Motif Instances": (df["LOST_MOTIF_INSTANCES"] > 0).sum(),
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    qgrs_change.plot(kind="bar", ax=axes[0], color=["#31a354", "#e6550d", "#9e9ac8"])
    axes[0].set_title(f"{region_label}: QGRS Count Change")
    axes[0].set_ylabel("Number of SNP-mutant sequences")
    axes[0].tick_params(axis="x", rotation=30)

    gscore_change.plot(kind="bar", ax=axes[1], color=["#3182bd", "#de2d26", "#9e9ac8"])
    axes[1].set_title(f"{region_label}: Max G-score Change")
    axes[1].tick_params(axis="x", rotation=30)

    motif_flags.plot(kind="bar", ax=axes[2], color=["#756bb1", "#bcbddc", "#41ab5d", "#fb6a4a"])
    axes[2].set_title(f"{region_label}: Motif-Level Change Flags")
    axes[2].tick_params(axis="x", rotation=35)

    fig.suptitle(f"{region_label} Wild vs Mutant Overview", y=1.03)
    fig.tight_layout()
    fig.savefig(out_dir / "01_overview_panels.png", bbox_inches="tight")
    plt.close(fig)

    save_series(effect_counts, out_dir.parent / "tables" / "effect_class_counts.tsv", "effect_class", "count")
    save_series(qgrs_change, out_dir.parent / "tables" / "qgrs_change_counts.tsv", "change_type", "count")
    save_series(gscore_change, out_dir.parent / "tables" / "gscore_change_counts.tsv", "change_type", "count")


def make_effect_class_bar(df, region_label, out_dir):
    counts = df["EFFECT_CLASS"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(counts.index, counts.values, color="#4c78a8")
    ax.set_title(f"{region_label}: Effect Class Counts")
    ax.set_xlabel("Count")
    fig.tight_layout()
    fig.savefig(out_dir / "02_effect_class_counts.png", bbox_inches="tight")
    plt.close(fig)


def make_delta_gscore_hist(df, region_label, out_dir):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(df["DELTA_MAX_GSCORE"], bins=30, color="#3182bd", edgecolor="white")
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(f"{region_label}: Distribution of Delta Max G-score")
    ax.set_xlabel("Mutant max g-score - Wild max g-score")
    ax.set_ylabel("Number of SNP-mutant sequences")
    fig.tight_layout()
    fig.savefig(out_dir / "03_delta_max_gscore_histogram.png", bbox_inches="tight")
    plt.close(fig)


def make_delta_qgrs_bar(df, region_label, out_dir):
    counts = df["DELTA_QGRS_COUNT"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(counts.index.astype(str), counts.values, color="#31a354")
    ax.set_title(f"{region_label}: Distribution of Delta QGRS Count")
    ax.set_xlabel("Mutant QGRS count - Wild QGRS count")
    ax.set_ylabel("Number of SNP-mutant sequences")
    fig.tight_layout()
    fig.savefig(out_dir / "04_delta_qgrs_count_distribution.png", bbox_inches="tight")
    plt.close(fig)


def make_gscore_scatter(df, region_label, out_dir):
    color_map = {"from_G": "#de2d26", "to_G": "#31a354"}
    fig, ax = plt.subplots(figsize=(7, 7))
    for cls, sub in df.groupby("G_CHANGE_CLASS"):
        ax.scatter(
            sub["REFERENCE_MAX_GSCORE"],
            sub["MUTANT_MAX_GSCORE"],
            alpha=0.65,
            s=26,
            label=cls,
            color=color_map.get(cls, "#636363"),
        )
    max_val = max(df["REFERENCE_MAX_GSCORE"].max(), df["MUTANT_MAX_GSCORE"].max())
    ax.plot([0, max_val], [0, max_val], linestyle="--", color="black", linewidth=1)
    ax.set_title(f"{region_label}: Wild vs Mutant Max G-score")
    ax.set_xlabel("Wild max g-score")
    ax.set_ylabel("Mutant max g-score")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "05_wild_vs_mutant_max_gscore_scatter.png", bbox_inches="tight")
    plt.close(fig)


def make_qgrs_scatter(df, region_label, out_dir):
    color_map = {"from_G": "#de2d26", "to_G": "#31a354"}
    fig, ax = plt.subplots(figsize=(7, 7))
    for cls, sub in df.groupby("G_CHANGE_CLASS"):
        ax.scatter(
            sub["REFERENCE_QGRS_COUNT"],
            sub["MUTANT_QGRS_COUNT"],
            alpha=0.65,
            s=26,
            label=cls,
            color=color_map.get(cls, "#636363"),
        )
    max_val = max(df["REFERENCE_QGRS_COUNT"].max(), df["MUTANT_QGRS_COUNT"].max())
    ax.plot([0, max_val], [0, max_val], linestyle="--", color="black", linewidth=1)
    ax.set_title(f"{region_label}: Wild vs Mutant QGRS Count")
    ax.set_xlabel("Wild QGRS count")
    ax.set_ylabel("Mutant QGRS count")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "06_wild_vs_mutant_qgrs_count_scatter.png", bbox_inches="tight")
    plt.close(fig)


def make_gchange_boxplot(df, region_label, out_dir):
    groups = []
    labels = []
    for label in ["from_G", "to_G"]:
        sub = df[df["G_CHANGE_CLASS"] == label]["DELTA_MAX_GSCORE"]
        if len(sub) > 0:
            groups.append(sub.values)
            labels.append(label)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(groups, tick_labels=labels, patch_artist=True, boxprops=dict(facecolor="#9ecae1"))
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(f"{region_label}: Delta Max G-score by G-change Class")
    ax.set_ylabel("Mutant max g-score - Wild max g-score")
    fig.tight_layout()
    fig.savefig(out_dir / "07_delta_gscore_by_g_change_class.png", bbox_inches="tight")
    plt.close(fig)


def make_top_gene_gscore(df, region_label, out_dir):
    gene_scores = (
        df.groupby("GENE")
        .agg(
            total_abs_delta_max_gscore=("DELTA_MAX_GSCORE", lambda s: s.abs().sum()),
            mean_delta_max_gscore=("DELTA_MAX_GSCORE", "mean"),
            snp_count=("GENE", "size"),
            reference_qgrs_wildtype=("REFERENCE_QGRS_COUNT", "max"),
        )
        .sort_values("total_abs_delta_max_gscore", ascending=False)
        .head(15)
    )
    save_table(gene_scores.reset_index(), out_dir.parent / "tables" / "top_genes_by_abs_gscore_change.tsv")

    plot_df = gene_scores.sort_values("total_abs_delta_max_gscore")
    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.barh(plot_df.index, plot_df["total_abs_delta_max_gscore"], color="#756bb1")
    qgrs_labels = [f"n={int(v)}" for v in plot_df["reference_qgrs_wildtype"]]
    ax.bar_label(bars, labels=qgrs_labels, padding=4, fontsize=9, color="#333333")
    ax.set_title(f"{region_label}: Top Genes by Total Absolute G-score Change")
    ax.set_xlabel("Sum of absolute delta max g-score")
    ax.set_ylabel("")
    key = (
        "Keyword — n: total QGRS motifs in the\n"
        "wild-type (reference) sequence for that gene in this region."
    )
    ax.text(
        0.99,
        0.01,
        key,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        linespacing=1.25,
        bbox=dict(
            boxstyle="round,pad=0.45",
            facecolor="#f8f8f0",
            edgecolor="#6a51a3",
            linewidth=1,
        ),
    )
    fig.tight_layout()
    fig.savefig(out_dir / "08_top_genes_by_abs_gscore_change.png", bbox_inches="tight")
    plt.close(fig)


def make_top_gene_qgrs(df, region_label, out_dir):
    gene_qgrs = (
        df.groupby("GENE")
        .agg(
            total_delta_qgrs_count=("DELTA_QGRS_COUNT", "sum"),
            abs_total_delta_qgrs=("DELTA_QGRS_COUNT", lambda s: s.abs().sum()),
            snp_count=("GENE", "size"),
        )
        .sort_values("abs_total_delta_qgrs", ascending=False)
        .head(15)
    )
    save_table(gene_qgrs.reset_index(), out_dir.parent / "tables" / "top_genes_by_qgrs_count_change.tsv")

    plot_df = gene_qgrs.sort_values("total_delta_qgrs_count")
    colors = ["#de2d26" if x < 0 else "#31a354" for x in plot_df["total_delta_qgrs_count"]]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(plot_df.index, plot_df["total_delta_qgrs_count"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_title(f"{region_label}: Top Genes by Net QGRS Count Change")
    ax.set_xlabel("Net delta QGRS count across SNPs")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(out_dir / "09_top_genes_by_qgrs_count_change.png", bbox_inches="tight")
    plt.close(fig)


def make_clnsig_effect(df, region_label, out_dir):
    effect_group = df["EFFECT_CLASS"].replace(
        {
            "higher_max_gscore": "higher_max_gscore",
            "lower_max_gscore": "lower_max_gscore",
            "no_summary_change": "no_summary_change",
        }
    )
    df = df.copy()
    df["EFFECT_GROUP"] = effect_group.where(
        effect_group.isin(["higher_max_gscore", "lower_max_gscore", "no_summary_change"]),
        "other_effects",
    )
    counts = (
        df.groupby(["CLNSIG_PRIMARY", "EFFECT_GROUP"]).size().reset_index(name="count")
    )
    top_categories = (
        counts.groupby("CLNSIG_PRIMARY")["count"].sum().sort_values(ascending=False).head(8).index
    )
    table = (
        counts[counts["CLNSIG_PRIMARY"].isin(top_categories)]
        .pivot(index="CLNSIG_PRIMARY", columns="EFFECT_GROUP", values="count")
        .fillna(0)
        .astype(int)
        .loc[top_categories]
    )
    table.to_csv(out_dir.parent / "tables" / "clinsig_by_effect_group.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(11, 6))
    table.plot(kind="bar", stacked=True, ax=ax, color=["#31a354", "#de2d26", "#9e9ac8", "#756bb1"])
    ax.set_title(f"{region_label}: ClinVar Categories by Effect Group")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "10_clinsig_by_effect_group.png", bbox_inches="tight")
    plt.close(fig)


def make_sorted_delta_line(df, region_label, out_dir):
    values = df["DELTA_MAX_GSCORE"].sort_values().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(range(1, len(values) + 1), values.values, color="#08519c", linewidth=2)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(f"{region_label}: Sorted Delta Max G-score Profile")
    ax.set_xlabel("Mutant SNP sequences ranked by delta max g-score")
    ax.set_ylabel("Delta max g-score")
    fig.tight_layout()
    fig.savefig(out_dir / "11_sorted_delta_max_gscore_profile.png", bbox_inches="tight")
    plt.close(fig)


def main():
    ensure_dirs()
    configure_plotting()

    for region_key, meta in REGIONS.items():
        df = load_region_df(region_key)
        region_dir = OUTPUT_ROOT / region_key
        fig_dir = region_dir / "figures"

        save_table(df, region_dir / "tables" / f"{region_key}_comparison_summary_copy.tsv")

        make_overview_figure(df, region_key, meta["label"], fig_dir)
        make_effect_class_bar(df, meta["label"], fig_dir)
        make_delta_gscore_hist(df, meta["label"], fig_dir)
        make_delta_qgrs_bar(df, meta["label"], fig_dir)
        make_gscore_scatter(df, meta["label"], fig_dir)
        make_qgrs_scatter(df, meta["label"], fig_dir)
        make_gchange_boxplot(df, meta["label"], fig_dir)
        make_top_gene_gscore(df, meta["label"], fig_dir)
        make_top_gene_qgrs(df, meta["label"], fig_dir)
        make_clnsig_effect(df, meta["label"], fig_dir)
        make_sorted_delta_line(df, meta["label"], fig_dir)

        print(f"Finished figure set for {region_key}")


if __name__ == "__main__":
    main()
