#!/usr/bin/env python3
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/Users/ritika/Desktop/pilot/results_cancer_context_figures/.matplotlib")

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path("/Users/ritika/Desktop/pilot")
COSMIC_PATH = PROJECT_ROOT / "gene" / "Cosmic_CancerGeneCensus_v103_GRCh38.tsv"
COMPARE_ROOT = PROJECT_ROOT / "compare_mutant_wild"
OUTPUT_ROOT = PROJECT_ROOT / "results_cancer_context_figures"
MPL_DIR = OUTPUT_ROOT / ".matplotlib"

REGIONS = {
    "tss": {
        "label": "TSS",
        "comparison": COMPARE_ROOT / "tss" / "tss_upstream500.comparison_summary.tsv",
    },
    "5utr": {
        "label": "5'UTR",
        "comparison": COMPARE_ROOT / "5utr" / "5utr.comparison_summary.tsv",
    },
    "3utr": {
        "label": "3'UTR",
        "comparison": COMPARE_ROOT / "3utr" / "3utr.comparison_summary.tsv",
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


def expand_multivalue(df, column):
    series = (
        df[column]
        .fillna("")
        .astype(str)
        .str.split(",")
        .explode()
        .str.strip()
        .replace("", "not_reported")
    )
    out = df.loc[series.index].copy()
    out[column] = series.values
    return out


def somatic_germline_group(row):
    somatic = str(row["SOMATIC"]).strip().lower()
    germline = str(row["GERMLINE"]).strip().lower()
    if somatic == "y" and germline == "y":
        return "somatic_and_germline"
    if somatic == "y" and germline != "y":
        return "somatic_only"
    if somatic != "y" and germline == "y":
        return "germline_only"
    return "neither_or_not_reported"


def load_region_df(region_key, cosmic_df):
    df = pd.read_csv(REGIONS[region_key]["comparison"], sep="\t")
    df["DELTA_MAX_GSCORE"] = pd.to_numeric(df["DELTA_MAX_GSCORE"])
    df["DELTA_QGRS_COUNT"] = pd.to_numeric(df["DELTA_QGRS_COUNT"])
    merged = df.merge(
        cosmic_df[
            [
                "GENE_SYMBOL",
                "NAME",
                "COSMIC_GENE_ID",
                "SOMATIC",
                "GERMLINE",
                "TISSUE_TYPE",
                "ROLE_IN_CANCER",
                "TUMOUR_TYPES_SOMATIC",
                "TUMOUR_TYPES_GERMLINE",
                "MOLECULAR_GENETICS",
                "TIER",
            ]
        ],
        left_on="GENE",
        right_on="GENE_SYMBOL",
        how="left",
    )
    merged["SOMATIC_GERMLINE_GROUP"] = merged.apply(somatic_germline_group, axis=1)
    return merged


def save_table(df, path):
    df.to_csv(path, sep="\t", index=False)


def make_role_count_plot(df, region_label, out_dir):
    role_df = expand_multivalue(df, "ROLE_IN_CANCER")
    counts = role_df["ROLE_IN_CANCER"].value_counts().head(12).sort_values()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(counts.index, counts.values, color="#4c78a8")
    ax.set_title(f"{region_label}: Variant Counts by Role in Cancer")
    ax.set_xlabel("Number of variant comparisons")
    fig.tight_layout()
    fig.savefig(out_dir / "01_role_in_cancer_variant_counts.png", bbox_inches="tight")
    plt.close(fig)
    save_table(counts.rename("count").reset_index().rename(columns={"index": "ROLE_IN_CANCER"}), out_dir.parent / "tables" / "role_in_cancer_variant_counts.tsv")


def make_role_delta_plot(df, region_label, out_dir):
    role_df = expand_multivalue(df, "ROLE_IN_CANCER")
    summary = (
        role_df.groupby("ROLE_IN_CANCER")
        .agg(mean_delta_max_gscore=("DELTA_MAX_GSCORE", "mean"), count=("ROLE_IN_CANCER", "size"))
        .query("count >= 5")
        .sort_values("mean_delta_max_gscore")
        .tail(12)
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["#de2d26" if x < 0 else "#31a354" for x in summary["mean_delta_max_gscore"]]
    ax.barh(summary.index, summary["mean_delta_max_gscore"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_title(f"{region_label}: Mean Delta Max G-score by Role in Cancer")
    ax.set_xlabel("Mean mutant - wild max g-score")
    fig.tight_layout()
    fig.savefig(out_dir / "02_role_in_cancer_mean_delta_gscore.png", bbox_inches="tight")
    plt.close(fig)
    save_table(summary.reset_index(), out_dir.parent / "tables" / "role_in_cancer_mean_delta_gscore.tsv")


def make_tissue_count_plot(df, region_label, out_dir):
    tissue_df = expand_multivalue(df, "TISSUE_TYPE")
    counts = tissue_df["TISSUE_TYPE"].value_counts().head(12).sort_values()
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(counts.index, counts.values, color="#756bb1")
    ax.set_title(f"{region_label}: Variant Counts by COSMIC Tissue Type")
    ax.set_xlabel("Number of variant comparisons")
    fig.tight_layout()
    fig.savefig(out_dir / "03_tissue_type_variant_counts.png", bbox_inches="tight")
    plt.close(fig)
    save_table(counts.rename("count").reset_index().rename(columns={"index": "TISSUE_TYPE"}), out_dir.parent / "tables" / "tissue_type_variant_counts.tsv")


def make_tissue_delta_plot(df, region_label, out_dir):
    tissue_df = expand_multivalue(df, "TISSUE_TYPE")
    summary = (
        tissue_df.groupby("TISSUE_TYPE")
        .agg(mean_delta_max_gscore=("DELTA_MAX_GSCORE", "mean"), count=("TISSUE_TYPE", "size"))
        .query("count >= 5")
        .sort_values("mean_delta_max_gscore")
        .tail(12)
    )
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#de2d26" if x < 0 else "#31a354" for x in summary["mean_delta_max_gscore"]]
    ax.barh(summary.index, summary["mean_delta_max_gscore"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_title(f"{region_label}: Mean Delta Max G-score by Tissue Type")
    ax.set_xlabel("Mean mutant - wild max g-score")
    fig.tight_layout()
    fig.savefig(out_dir / "04_tissue_type_mean_delta_gscore.png", bbox_inches="tight")
    plt.close(fig)
    save_table(summary.reset_index(), out_dir.parent / "tables" / "tissue_type_mean_delta_gscore.tsv")


def make_somatic_germline_counts(df, region_label, out_dir):
    counts = df["SOMATIC_GERMLINE_GROUP"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(counts.index, counts.values, color="#f28e2b")
    ax.set_title(f"{region_label}: Variant Counts by Somatic/Germline Group")
    ax.set_xlabel("Number of variant comparisons")
    fig.tight_layout()
    fig.savefig(out_dir / "05_somatic_germline_variant_counts.png", bbox_inches="tight")
    plt.close(fig)
    save_table(counts.rename("count").reset_index().rename(columns={"index": "SOMATIC_GERMLINE_GROUP"}), out_dir.parent / "tables" / "somatic_germline_variant_counts.tsv")


def make_somatic_germline_delta(df, region_label, out_dir):
    summary = (
        df.groupby("SOMATIC_GERMLINE_GROUP")
        .agg(mean_delta_max_gscore=("DELTA_MAX_GSCORE", "mean"), count=("SOMATIC_GERMLINE_GROUP", "size"))
        .sort_values("mean_delta_max_gscore")
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#de2d26" if x < 0 else "#31a354" for x in summary["mean_delta_max_gscore"]]
    ax.barh(summary.index, summary["mean_delta_max_gscore"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_title(f"{region_label}: Mean Delta Max G-score by Somatic/Germline Group")
    ax.set_xlabel("Mean mutant - wild max g-score")
    fig.tight_layout()
    fig.savefig(out_dir / "06_somatic_germline_mean_delta_gscore.png", bbox_inches="tight")
    plt.close(fig)
    save_table(summary.reset_index(), out_dir.parent / "tables" / "somatic_germline_mean_delta_gscore.tsv")


def make_tier_counts(df, region_label, out_dir):
    counts = df["TIER"].fillna("not_reported").astype(str).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(counts.index, counts.values, color="#59a14f")
    ax.set_title(f"{region_label}: Variant Counts by COSMIC Tier")
    ax.set_xlabel("COSMIC tier")
    ax.set_ylabel("Number of variant comparisons")
    fig.tight_layout()
    fig.savefig(out_dir / "07_cosmic_tier_counts.png", bbox_inches="tight")
    plt.close(fig)
    save_table(counts.rename("count").reset_index().rename(columns={"index": "TIER"}), out_dir.parent / "tables" / "cosmic_tier_counts.tsv")


def make_molecular_genetics_counts(df, region_label, out_dir):
    mg_df = expand_multivalue(df, "MOLECULAR_GENETICS")
    counts = mg_df["MOLECULAR_GENETICS"].value_counts().head(10).sort_values()
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(counts.index, counts.values, color="#e15759")
    ax.set_title(f"{region_label}: Variant Counts by Molecular Genetics")
    ax.set_xlabel("Number of variant comparisons")
    fig.tight_layout()
    fig.savefig(out_dir / "08_molecular_genetics_counts.png", bbox_inches="tight")
    plt.close(fig)
    save_table(counts.rename("count").reset_index().rename(columns={"index": "MOLECULAR_GENETICS"}), out_dir.parent / "tables" / "molecular_genetics_counts.tsv")


def build_text_summary(region_key, region_label, df):
    role_df = expand_multivalue(df, "ROLE_IN_CANCER")
    tissue_df = expand_multivalue(df, "TISSUE_TYPE")
    role_counts = role_df["ROLE_IN_CANCER"].value_counts()
    tissue_counts = tissue_df["TISSUE_TYPE"].value_counts()
    som_counts = df["SOMATIC_GERMLINE_GROUP"].value_counts()
    role_delta = (
        role_df.groupby("ROLE_IN_CANCER")["DELTA_MAX_GSCORE"].mean().sort_values()
    )
    tissue_delta = (
        tissue_df.groupby("TISSUE_TYPE")["DELTA_MAX_GSCORE"].mean().sort_values()
    )
    lines = []
    lines.append(f"{region_label} COSMIC Cancer Context Summary")
    lines.append("")
    lines.append(f"Total variant comparisons in this region: {len(df)}")
    lines.append(f"Unique genes represented: {df['GENE'].nunique()}")
    lines.append(f"Mean delta max g-score: {round(df['DELTA_MAX_GSCORE'].mean(), 3)}")
    lines.append("")
    lines.append("Most frequent role-in-cancer categories")
    for role, count in role_counts.head(8).items():
        lines.append(f"- {role}: {count}")
    lines.append("")
    lines.append("Role-in-cancer categories with most negative mean delta max g-score")
    for role, value in role_delta.head(5).items():
        lines.append(f"- {role}: {round(value, 3)}")
    lines.append("")
    lines.append("Most frequent tissue-type categories")
    for tissue, count in tissue_counts.head(8).items():
        lines.append(f"- {tissue}: {count}")
    lines.append("")
    lines.append("Tissue-type categories with most negative mean delta max g-score")
    for tissue, value in tissue_delta.head(5).items():
        lines.append(f"- {tissue}: {round(value, 3)}")
    lines.append("")
    lines.append("Somatic/germline group counts")
    for group, count in som_counts.items():
        lines.append(f"- {group}: {count}")
    lines.append("")
    lines.append("Interpretation")
    lines.append(f"- This summary joins the {region_label} wild-vs-mutant comparison results to COSMIC cancer-gene annotations.")
    lines.append("- Counts reflect variant comparisons, not just unique genes, so genes with multiple qualifying SNPs can contribute multiple rows.")
    lines.append("- Role-in-cancer and tissue-type fields were expanded when multiple labels were present in the COSMIC source file.")
    lines.append("- Somatic/germline grouping was derived directly from the COSMIC SOMATIC and GERMLINE columns.")
    return "\n".join(lines) + "\n"


def main():
    ensure_dirs()
    configure_plotting()

    cosmic_df = pd.read_csv(COSMIC_PATH, sep="\t")

    summary_lines = []
    summary_lines.append("COSMIC Cancer Context Figure Summary")
    summary_lines.append("")
    summary_lines.append("This output set joins the wild-vs-mutant comparison results to annotation fields from the COSMIC Cancer Gene Census, including tissue type, role in cancer, and somatic/germline status.")
    summary_lines.append("")

    for region_key, meta in REGIONS.items():
        region_dir = OUTPUT_ROOT / region_key
        fig_dir = region_dir / "figures"
        table_dir = region_dir / "tables"
        df = load_region_df(region_key, cosmic_df)

        save_table(df, table_dir / f"{region_key}_cosmic_annotated_comparison.tsv")

        make_role_count_plot(df, meta["label"], fig_dir)
        make_role_delta_plot(df, meta["label"], fig_dir)
        make_tissue_count_plot(df, meta["label"], fig_dir)
        make_tissue_delta_plot(df, meta["label"], fig_dir)
        make_somatic_germline_counts(df, meta["label"], fig_dir)
        make_somatic_germline_delta(df, meta["label"], fig_dir)
        make_tier_counts(df, meta["label"], fig_dir)
        make_molecular_genetics_counts(df, meta["label"], fig_dir)

        text = build_text_summary(region_key, meta["label"], df)
        (region_dir / f"{region_key}_cosmic_context_summary.txt").write_text(text)

        summary_lines.append(f"{meta['label']}: {len(df)} variant comparisons, {df['GENE'].nunique()} unique genes")

        print(f"Finished COSMIC context figure set for {region_key}")

    (OUTPUT_ROOT / "cosmic_context_results_summary.txt").write_text("\n".join(summary_lines) + "\n")


if __name__ == "__main__":
    main()
