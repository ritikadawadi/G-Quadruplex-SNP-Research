#!/usr/bin/env python3
import csv
from pathlib import Path


PROJECT_ROOT = Path("/Users/ritika/Desktop/pilot")
OUTPUT_ROOT = PROJECT_ROOT / "results_changes_distinct"

REGIONS = {
    "tss": {
        "label": "TSS",
        "comparison": PROJECT_ROOT / "compare_mutant_wild" / "tss" / "tss_upstream500.comparison_summary.tsv",
        "reference_fasta": PROJECT_ROOT / "data" / "fasta" / "tss_upstream500.reference.fa",
        "mutant_fasta": PROJECT_ROOT / "mutant_fasta" / "tss" / "tss_upstream500.mutant.fa",
    },
    "5utr": {
        "label": "5'UTR",
        "comparison": PROJECT_ROOT / "compare_mutant_wild" / "5utr" / "5utr.comparison_summary.tsv",
        "reference_fasta": PROJECT_ROOT / "data" / "fasta" / "5utr.reference.fa",
        "mutant_fasta": PROJECT_ROOT / "mutant_fasta" / "5utr" / "5utr.mutant.fa",
    },
    "3utr": {
        "label": "3'UTR",
        "comparison": PROJECT_ROOT / "compare_mutant_wild" / "3utr" / "3utr.comparison_summary.tsv",
        "reference_fasta": PROJECT_ROOT / "data" / "fasta" / "3utr.reference.fa",
        "mutant_fasta": PROJECT_ROOT / "mutant_fasta" / "3utr" / "3utr.mutant.fa",
    },
}


def read_fasta(path):
    records = {}
    header = None
    seq_parts = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records[header] = "".join(seq_parts)
    return records


def load_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def context_with_highlight(sequence, pos1, context=35):
    idx0 = int(pos1) - 1
    start = max(0, idx0 - context)
    end = min(len(sequence), idx0 + context + 1)
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(sequence) else ""
    return prefix + sequence[start:idx0] + "[" + sequence[idx0] + "]" + sequence[idx0 + 1 : end] + suffix


def format_entry(row, ref_seq, mut_seq):
    wt_len = len(row["REFERENCE_TOP_SEQUENCE"])
    mut_len = len(row["MUTANT_TOP_SEQUENCE"])
    reasons = []
    if int(row["DELTA_MAX_GSCORE"]) > 0:
        reasons.append("g-score increase")
    if int(row["DELTA_MAX_GSCORE"]) < 0:
        reasons.append("g-score decrease")
    if mut_len > wt_len:
        reasons.append("QGRS length increase")
    if mut_len < wt_len:
        reasons.append("QGRS length decrease")
    reason_text = ", ".join(reasons) if reasons else "selected change"

    clinvar_id = row["ALLELEID"] or row["VARIANT_ID"]
    variant_desc = f"{row['VARIANT_CHROM']}:{row['VARIANT_POS_1BASED']}"
    change_desc = f"{row['REF_ALLELE_GENOMIC']}>{row['ALT_ALLELE_GENOMIC']}"

    wt_context = context_with_highlight(ref_seq, row["TRANSCRIPT_POS1"])
    mut_context = context_with_highlight(mut_seq, row["TRANSCRIPT_POS1"])

    lines = []
    lines.append(f"{row['GENE']} | {change_desc} | ClinVar:{clinvar_id} | chr{variant_desc}")
    lines.append(
        f"WT  G-score: {row['REFERENCE_MAX_GSCORE']} | QGRS length: {wt_len} | Sequence: {wt_context}"
    )
    lines.append(
        f"MUT G-score: {row['MUTANT_MAX_GSCORE']} | QGRS length: {mut_len} | Sequence: {mut_context}"
    )
    lines.append(f"ΔG-score: {row['DELTA_MAX_GSCORE']}")
    lines.append(f"ΔQGRS length: {mut_len - wt_len}")
    lines.append(f"Reason included: {reason_text}")
    lines.append(f"G-change class: {row['G_CHANGE_CLASS']}")
    lines.append(f"WT top QGRS:  {row['REFERENCE_TOP_SEQUENCE']}")
    lines.append(f"MUT top QGRS: {row['MUTANT_TOP_SEQUENCE']}")
    lines.append("")
    return "\n".join(lines)


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    for region_key, config in REGIONS.items():
        region_dir = OUTPUT_ROOT / region_key
        region_dir.mkdir(parents=True, exist_ok=True)

        comparison_rows = load_rows(config["comparison"])
        reference_fasta = read_fasta(config["reference_fasta"])
        mutant_fasta = read_fasta(config["mutant_fasta"])

        selected = []
        for row in comparison_rows:
            wt_len = len(row["REFERENCE_TOP_SEQUENCE"])
            mut_len = len(row["MUTANT_TOP_SEQUENCE"])
            if int(row["DELTA_MAX_GSCORE"]) != 0 or mut_len != wt_len:
                selected.append(row)

        out_path = region_dir / f"{region_key}_gscore_or_length_change_report.txt"
        with out_path.open("w") as handle:
            handle.write("=" * 100 + "\n")
            handle.write(f"QGRS DISTINCT CHANGE REPORT - {config['label']}\n")
            handle.write("=" * 100 + "\n\n")
            handle.write(f"Total qualifying variants: {len(selected)}\n")
            handle.write(f"Unique genes: {len({row['GENE'] for row in selected})}\n")
            handle.write(
                f"Rule: include cases where mutant max g-score differs from wild max g-score OR mutant top QGRS length differs from wild top QGRS length.\n\n"
            )

            for row in selected:
                ref_header = row["REFERENCE_HEADER"]
                mut_header = row["MUTANT_HEADER"]
                ref_seq = reference_fasta[ref_header]
                mut_seq = mutant_fasta[mut_header]
                handle.write(format_entry(row, ref_seq, mut_seq))
                handle.write("-" * 100 + "\n\n")

        print(f"Wrote {out_path}")
        print(f"Qualifying variants for {region_key}: {len(selected)}")


if __name__ == "__main__":
    main()
