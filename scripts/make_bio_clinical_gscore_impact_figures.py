#!/usr/bin/env python3
"""
ClinVar clinical significance + COSMIC cancer role vs QGRS max g-score change.

Poster outputs (300 DPI) go to results_compare_figures/{tss,5utr,3utr}/figures/
  13_clinsig_pct_gscore_changed.png
  14_clinsig_mean_abs_delta_when_changed.png
  15_clinsig_stacked_pct_gscore_change.png
  16_cosmic_role_pct_gscore_changed.png
  17_cosmic_role_mean_abs_delta_when_changed.png
  18_clndn_pct_gscore_changed.png
  19_clndn_mean_abs_delta_when_changed.png
  20_clndn_stacked_pct_gscore_change.png

Tables: results_compare_figures/{region}/tables/bio_clinical_*.tsv
"""
from __future__ import annotations

import os
import re
import textwrap
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent.parent / "results_compare_figures" / ".matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_ROOT = PROJECT_ROOT / "results_cancer_context_figures"
COMPARE_FIG_ROOT = PROJECT_ROOT / "results_compare_figures"
POSTER_DPI = 300

REGIONS = {
    "tss": {
        "label": "TSS (upstream 500 bp)",
        "annotated": INPUT_ROOT / "tss" / "tables" / "tss_cosmic_annotated_comparison.tsv",
    },
    "5utr": {
        "label": "5′UTR",
        "annotated": INPUT_ROOT / "5utr" / "tables" / "5utr_cosmic_annotated_comparison.tsv",
    },
    "3utr": {
        "label": "3′UTR",
        "annotated": INPUT_ROOT / "3utr" / "tables" / "3utr_cosmic_annotated_comparison.tsv",
    },
}


def configure_poster_style():
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": POSTER_DPI,
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )


def expand_multivalue(df: pd.DataFrame, column: str) -> pd.DataFrame:
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


def simplify_clnsig(raw) -> str:
    if pd.isna(raw) or str(raw).strip() == "":
        return "Not reported"
    low = str(raw).strip().lower()
    if "conflicting" in low:
        return "Conflicting"
    if "uncertain" in low:
        return "Uncertain significance"
    if "likely_benign" in low or "likely benign" in low:
        return "Likely benign"
    if "benign" in low and "pathogenic" not in low:
        return "Benign"
    if "likely_pathogenic" in low or "likely pathogenic" in low:
        return "Likely pathogenic"
    if "pathogenic" in low:
        return "Pathogenic"
    return "Other"


def role_short(x: str) -> str:
    x = str(x).strip()
    if not x or x == "not_reported":
        return "Not reported"
    return x[:48] + ("…" if len(x) > 48 else "")


JUNK_CLNDN = frozenset(
    {
        "not_specified",
        "not specified",
        "not_provided",
        "not provided",
        "-",
        "see_cases",
        "see cases",
        "all_submitters",
        "all submitters",
        "nan",
        "none",
        "",
    }
)


def clean_clndn_term(s: str) -> str:
    s = str(s).strip()
    if not s:
        return ""
    s = re.sub(r"_+", " ", s).strip()
    low = s.lower()
    if low in JUNK_CLNDN or "not_specified" in low or "not provided" in low:
        return ""
    if len(s) < 2:
        return ""
    return s


def expand_clndn(df: pd.DataFrame) -> pd.DataFrame:
    """One row per variant×disease term (CLNDN pipe-separated)."""
    if "CLNDN" not in df.columns:
        return pd.DataFrame()
    out_rows = []
    for _, row in df.iterrows():
        raw = row["CLNDN"]
        if pd.isna(raw) or str(raw).strip() == "":
            continue
        parts = re.split(r"\s*\|\s*", str(raw))
        seen: set[str] = set()
        for p in parts:
            term = clean_clndn_term(p)
            if not term:
                continue
            key = term.casefold()
            if key in seen:
                continue
            seen.add(key)
            r = row.copy()
            r["CLNDN_TERM"] = term
            out_rows.append(r)
    if not out_rows:
        return pd.DataFrame()
    return pd.DataFrame(out_rows)


def summarize_category(df: pd.DataFrame, col: str, min_n: int = 8) -> pd.DataFrame:
    work = df[[col, "DELTA_MAX_GSCORE", "abs_delta", "g_changed"]].copy()
    work = work[work[col].notna() & (work[col].astype(str).str.strip() != "")]
    rows = []
    for cat, g in work.groupby(col, dropna=False):
        n = len(g)
        if n < min_n:
            continue
        ch = int(g["g_changed"].sum())
        sub = g[g["g_changed"]]
        rows.append(
            {
                col: cat,
                "n_variants": n,
                "n_gscore_changed": ch,
                "pct_gscore_changed": 100.0 * ch / n if n else 0.0,
                "mean_abs_delta_all": float(g["abs_delta"].mean()),
                "mean_abs_delta_when_changed": float(sub["abs_delta"].mean())
                if len(sub)
                else float("nan"),
                "mean_signed_delta": float(g["DELTA_MAX_GSCORE"].mean()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pct_gscore_changed", ascending=False)


def _wrap_ytick_labels(labels: list, wrap_width: int | None) -> list:
    if not wrap_width:
        return labels
    return [textwrap.fill(str(lab), wrap_width) for lab in labels]


def plot_pct_changed(
    summary: pd.DataFrame,
    col: str,
    title: str,
    outfile: Path,
    *,
    y_label_wrap: int | None = None,
) -> None:
    d = summary.sort_values("pct_gscore_changed").tail(14)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(11, max(5.5, 0.4 * len(d) + 1.2)))
    y = np.arange(len(d))
    ax.barh(y, d["pct_gscore_changed"], color="#3182bd", height=0.72, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y, _wrap_ytick_labels(list(d[col]), y_label_wrap))
    ax.set_xlabel("% variants with Δ max g-score ≠ 0")
    ax.set_title(title)
    xmax = max(15.0, float(d["pct_gscore_changed"].max()) * 1.12)
    ax.set_xlim(0, min(100.0, xmax))
    ax.grid(axis="x", alpha=0.35)
    fig.tight_layout()
    fig.savefig(outfile, dpi=POSTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_mean_abs_when_changed(
    summary: pd.DataFrame,
    col: str,
    title: str,
    outfile: Path,
    *,
    y_label_wrap: int | None = None,
) -> None:
    d = summary.dropna(subset=["mean_abs_delta_when_changed"]).copy()
    d = d[d["n_gscore_changed"] >= 3]
    d = d.sort_values("mean_abs_delta_when_changed").tail(14)
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(11, max(5.5, 0.4 * len(d) + 1.2)))
    y = np.arange(len(d))
    ax.barh(y, d["mean_abs_delta_when_changed"], color="#e6550d", height=0.72, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y, _wrap_ytick_labels(list(d[col]), y_label_wrap))
    ax.set_xlabel("Mean |Δ max g-score| (among variants with change)")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.35)
    fig.tight_layout()
    fig.savefig(outfile, dpi=POSTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_grouped_clinsig_effect(df: pd.DataFrame, region_label: str, outfile: Path) -> None:
    vc = df.groupby("CLNSIG_SIMPLE")["g_changed"].agg(["mean", "count"])
    cats = vc[vc["count"] >= 5].index.tolist()
    if not cats:
        return
    cats = sorted(cats, key=lambda c: -vc.loc[c, "mean"])
    pct_changed = [100 * float(df[df["CLNSIG_SIMPLE"] == c]["g_changed"].mean()) for c in cats]
    pct_stable = [100 - p for p in pct_changed]

    fig, ax = plt.subplots(figsize=(12, 6.2))
    x = np.arange(len(cats))
    w = 0.62
    ax.bar(x, pct_changed, w, label="Δ max g-score ≠ 0", color="#31a354", edgecolor="white", linewidth=0.6)
    ax.bar(
        x,
        pct_stable,
        w,
        bottom=pct_changed,
        label="Δ max g-score = 0",
        color="#d9d9d9",
        edgecolor="white",
        linewidth=0.6,
    )
    ax.set_xticks(x, [c.replace(" ", "\n") for c in cats], fontsize=10)
    ax.set_ylabel("Percent of variants")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92)
    ax.set_title(f"{region_label}: G-score change by ClinVar clinical significance")
    ax.grid(axis="y", alpha=0.35)
    fig.tight_layout()
    fig.savefig(outfile, dpi=POSTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_grouped_disease_effect(
    df: pd.DataFrame,
    region_label: str,
    outfile: Path,
    *,
    top_n: int = 10,
    min_n: int = 8,
) -> None:
    vc = df.groupby("CLNDN_TERM")["g_changed"].agg(["mean", "count"])
    vc = vc[vc["count"] >= min_n]
    if vc.empty:
        return
    cats = vc.sort_values("count", ascending=False).head(top_n).index.tolist()
    cats = sorted(cats, key=lambda c: -float(vc.loc[c, "mean"]))
    pct_changed = [100 * float(df[df["CLNDN_TERM"] == c]["g_changed"].mean()) for c in cats]
    pct_stable = [100 - p for p in pct_changed]

    fig, ax = plt.subplots(figsize=(14, 6.2))
    x = np.arange(len(cats))
    w = 0.62
    ax.bar(x, pct_changed, w, label="Δ max g-score ≠ 0", color="#31a354", edgecolor="white", linewidth=0.6)
    ax.bar(
        x,
        pct_stable,
        w,
        bottom=pct_changed,
        label="Δ max g-score = 0",
        color="#d9d9d9",
        edgecolor="white",
        linewidth=0.6,
    )
    xlabs = [textwrap.fill(str(c), 22) for c in cats]
    ax.set_xticks(x, xlabs, fontsize=9)
    ax.set_ylabel("Percent of variants")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92)
    ax.set_title(f"{region_label}: G-score change by ClinVar disease (CLNDN)")
    ax.grid(axis="y", alpha=0.35)
    fig.tight_layout()
    fig.savefig(outfile, dpi=POSTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_interpretation_note(region_label: str, out_dir: Path) -> None:
    path = out_dir / "bio_clinical_gscore_README.txt"
    body = (
        f"{region_label} — figures 13–20 (bio/clinical vs G-score)\n\n"
        "• 13 / 16 / 18: % of variants in each bucket with any Δ max g-score.\n"
        "• 14 / 17 / 19: mean |mutant−wild| among variants that do change.\n"
        "• 15 / 20: stacked bars (green = changed; grey = unchanged).\n\n"
        "CLNSIG, CLNDN (disease), and COSMIC roles annotate variants/genes, not G4 structure directly.\n"
        "COSMIC roles are comma-exploded; CLNDN is pipe-exploded — one variant may appear in multiple bins.\n"
        "Association only.\n"
    )
    path.write_text(body, encoding="utf-8")


def process_region(region_key: str) -> None:
    meta = REGIONS[region_key]
    df = pd.read_csv(meta["annotated"], sep="\t")
    df["DELTA_MAX_GSCORE"] = pd.to_numeric(df["DELTA_MAX_GSCORE"], errors="coerce").fillna(0.0)
    df["abs_delta"] = df["DELTA_MAX_GSCORE"].abs()
    df["g_changed"] = df["abs_delta"] > 0
    df["CLNSIG_SIMPLE"] = df["CLNSIG"].apply(simplify_clnsig) if "CLNSIG" in df.columns else "Not reported"

    fig_dir = COMPARE_FIG_ROOT / region_key / "figures"
    tbl_dir = COMPARE_FIG_ROOT / region_key / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    cl_summary = summarize_category(df, "CLNSIG_SIMPLE", min_n=5)
    snap_rows = []

    if not cl_summary.empty:
        cl_summary.to_csv(tbl_dir / "bio_clinical_clinsig_gscore_impact.tsv", sep="\t", index=False)
        plot_pct_changed(
            cl_summary,
            "CLNSIG_SIMPLE",
            f'{meta["label"]}: % variants with G-score change (ClinVar)',
            fig_dir / "13_clinsig_pct_gscore_changed.png",
        )
        plot_mean_abs_when_changed(
            cl_summary,
            "CLNSIG_SIMPLE",
            f'{meta["label"]}: |Δ max g-score| when changed (ClinVar)',
            fig_dir / "14_clinsig_mean_abs_delta_when_changed.png",
        )

    plot_grouped_clinsig_effect(
        df, meta["label"], fig_dir / "15_clinsig_stacked_pct_gscore_change.png"
    )

    role_summary = pd.DataFrame()
    if "ROLE_IN_CANCER" in df.columns:
        role_df = expand_multivalue(df, "ROLE_IN_CANCER")
        role_df["ROLE_SHORT"] = role_df["ROLE_IN_CANCER"].map(role_short)
        role_summary = summarize_category(role_df, "ROLE_SHORT", min_n=8)
        if not role_summary.empty:
            role_summary.to_csv(
                tbl_dir / "bio_clinical_cosmic_role_gscore_impact.tsv", sep="\t", index=False
            )
            plot_pct_changed(
                role_summary,
                "ROLE_SHORT",
                f'{meta["label"]}: % variants with G-score change (COSMIC role)',
                fig_dir / "16_cosmic_role_pct_gscore_changed.png",
            )
            plot_mean_abs_when_changed(
                role_summary,
                "ROLE_SHORT",
                f'{meta["label"]}: |Δ max g-score| when changed (COSMIC role)',
                fig_dir / "17_cosmic_role_mean_abs_delta_when_changed.png",
            )

    clnd_summary = pd.DataFrame()
    if "CLNDN" in df.columns:
        clnd_df = expand_clndn(df)
        if not clnd_df.empty:
            clnd_summary = summarize_category(clnd_df, "CLNDN_TERM", min_n=8)
            if not clnd_summary.empty:
                clnd_summary.to_csv(
                    tbl_dir / "bio_clinical_clndn_gscore_impact.tsv", sep="\t", index=False
                )
                plot_pct_changed(
                    clnd_summary,
                    "CLNDN_TERM",
                    f'{meta["label"]}: % variants with G-score change (ClinVar CLNDN)',
                    fig_dir / "18_clndn_pct_gscore_changed.png",
                    y_label_wrap=38,
                )
                plot_mean_abs_when_changed(
                    clnd_summary,
                    "CLNDN_TERM",
                    f'{meta["label"]}: |Δ max g-score| when changed (ClinVar CLNDN)',
                    fig_dir / "19_clndn_mean_abs_delta_when_changed.png",
                    y_label_wrap=38,
                )
            plot_grouped_disease_effect(
                clnd_df,
                meta["label"],
                fig_dir / "20_clndn_stacked_pct_gscore_change.png",
                top_n=10,
                min_n=8,
            )

    for _, r in cl_summary.iterrows():
        snap_rows.append({"region": region_key, "axis": "ClinVar_CLNSIG", **r.to_dict()})
    for _, r in role_summary.iterrows():
        snap_rows.append({"region": region_key, "axis": "COSMIC_role", **r.to_dict()})
    for _, r in clnd_summary.iterrows():
        snap_rows.append({"region": region_key, "axis": "ClinVar_CLNDN", **r.to_dict()})
    if snap_rows:
        pd.DataFrame(snap_rows).to_csv(
            tbl_dir / "bio_clinical_combined_summary.tsv", sep="\t", index=False
        )

    write_interpretation_note(meta["label"], tbl_dir)
    print(f"Wrote poster figures under {fig_dir}")


def main():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")
    configure_poster_style()
    COMPARE_FIG_ROOT.mkdir(parents=True, exist_ok=True)
    for key in REGIONS:
        apath = REGIONS[key]["annotated"]
        if not apath.is_file():
            raise FileNotFoundError(apath)
        process_region(key)


if __name__ == "__main__":
    main()
