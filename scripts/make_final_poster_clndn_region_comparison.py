#!/usr/bin/env python3
from pathlib import Path
import os
import textwrap

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "results_compare_figures" / "final poster figures"

JUNK = {"not provided", "not specified", "see cases", "not_provided", "not_specified", "see_cases"}

REGIONS = {
    "tss": {
        "label": "Upstream Region",
        "table": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "bio_clinical_clndn_gscore_impact.tsv",
        "summary": PROJECT_ROOT / "results_compare_figures" / "tss" / "tables" / "tss_comparison_summary_copy.tsv",
        "color": "#e6550d",
    },
    "5utr": {
        "label": "5\u2032UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "bio_clinical_clndn_gscore_impact.tsv",
        "summary": PROJECT_ROOT / "results_compare_figures" / "5utr" / "tables" / "5utr_comparison_summary_copy.tsv",
        "color": "#2ca25f",
    },
    "3utr": {
        "label": "3\u2032UTR",
        "table": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "bio_clinical_clndn_gscore_impact.tsv",
        "summary": PROJECT_ROOT / "results_compare_figures" / "3utr" / "tables" / "3utr_comparison_summary_copy.tsv",
        "color": "#756bb1",
    },
}


def compute_from_summary(summary_path: Path, top_n: int = 5) -> pd.DataFrame:
    """Build CLNDN stats directly from the comparison summary."""
    df = pd.read_csv(summary_path, sep="\t")
    rows = []
    for _, r in df.iterrows():
        clndn = str(r.get("CLNDN", ""))
        delta = r["DELTA_MAX_GSCORE"]
        for term in clndn.split("|"):
            term = term.strip().replace("_", " ")
            if term and term.lower() not in JUNK:
                rows.append({"term": term, "delta": delta})
    if not rows:
        return pd.DataFrame(columns=["CLNDN_TERM", "pct_gscore_changed", "n_variants", "label"])
    ex = pd.DataFrame(rows)
    g = ex.groupby("term").agg(
        n_variants=("delta", "count"),
        n_changed=("delta", lambda s: (s != 0).sum()),
    ).reset_index()
    g = g[g["n_variants"] >= 3]
    g["pct_gscore_changed"] = 100.0 * g["n_changed"] / g["n_variants"]
    g = g.sort_values(["pct_gscore_changed", "n_variants"], ascending=[False, False]).head(top_n)
    g = g.rename(columns={"term": "CLNDN_TERM"})
    g["label"] = g["CLNDN_TERM"].map(lambda x: textwrap.fill(str(x), width=24))
    return g


def load_top_terms(meta: dict, top_n: int = 5) -> pd.DataFrame:
    df = pd.read_csv(meta["table"], sep="\t")
    if len(df) >= top_n:
        df = df.sort_values(["pct_gscore_changed", "n_variants"], ascending=[False, False]).head(top_n).copy()
        df["label"] = df["CLNDN_TERM"].map(lambda x: textwrap.fill(str(x), width=24))
        return df
    return compute_from_summary(meta["summary"], top_n)


def add_panel(ax, df: pd.DataFrame, title: str, color: str) -> None:
    x = range(len(df))
    ax.bar(x, df["pct_gscore_changed"], color=color, edgecolor="white", linewidth=0.8, width=0.65)
    ax.set_xticks(list(x), df["label"], rotation=45, ha="right", fontsize=13)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of variants with\nmax G-score change", fontsize=16)
    ax.set_title(title, fontweight="bold", fontsize=20)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", alpha=0.25)
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            i,
            min(float(row["pct_gscore_changed"]) + 2, 97),
            f'n={int(row["n_variants"])}',
            ha="center",
            va="bottom",
            fontsize=12,
            color="#333333",
        )
    ax.text(0.97, 0.97, "n = no. of 'G'-related\nSNVs per disease",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color="#888888",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)

    for ax, meta in zip(axes, REGIONS.values()):
        df = load_top_terms(meta, top_n=5)
        add_panel(ax, df, meta["label"], meta["color"])

    fig.suptitle(
        "ClinVar disease terms most often linked to max G-score change",
        fontsize=24,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.01,
        "Top disease terms per region ranked by the percent of 'G'-related variants whose wild-type vs mutant sequence changed the strongest predicted G-score.",
        ha="center",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    outfile = OUTPUT_DIR / "combined_clndn_pct_gscore_changed_5utr_3utr.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
