#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def load_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_int(value):
    return int(value) if value not in (None, "", "NA") else 0


def motif_signature(row):
    return (
        row["QGRS_SEQUENCE"],
        safe_int(row["GSCORE"]),
        safe_int(row["TETRADS"]),
        safe_int(row["LENGTH"]),
    )


def best_hit(hits):
    if not hits:
        return None
    return max(
        hits,
        key=lambda row: (
            safe_int(row["GSCORE"]),
            safe_int(row["TETRADS"]),
            safe_int(row["LENGTH"]),
            row["QGRS_SEQUENCE"],
        ),
    )


def build_motif_delta_rows(
    manifest_row,
    ref_hits,
    mut_hits,
    reference_header,
    mutant_header,
):
    ref_counter = Counter(motif_signature(row) for row in ref_hits)
    mut_counter = Counter(motif_signature(row) for row in mut_hits)

    delta_rows = []
    for signature, count in sorted((mut_counter - ref_counter).items()):
        sequence, gscore, tetrads, length = signature
        delta_rows.append(
            {
                "MUTANT_HEADER": mutant_header,
                "REFERENCE_HEADER": reference_header,
                "GENE": manifest_row["GENE"],
                "TRANSCRIPT": manifest_row["TRANSCRIPT"],
                "STRAND": manifest_row["STRAND"],
                "REGION": manifest_row["REGION"],
                "VARIANT_CHROM": manifest_row["VARIANT_CHROM"],
                "VARIANT_POS_1BASED": manifest_row["VARIANT_POS_1BASED"],
                "REF_ALLELE_GENOMIC": manifest_row["REF_ALLELE_GENOMIC"],
                "ALT_ALLELE_GENOMIC": manifest_row["ALT_ALLELE_GENOMIC"],
                "REF_ALLELE_TRANSCRIPT": manifest_row["REF_ALLELE_TRANSCRIPT"],
                "ALT_ALLELE_TRANSCRIPT": manifest_row["ALT_ALLELE_TRANSCRIPT"],
                "TRANSCRIPT_POS1": manifest_row["TRANSCRIPT_POS1"],
                "G_CHANGE_CLASS": manifest_row["G_CHANGE_CLASS"],
                "DELTA_TYPE": "gained",
                "DELTA_COUNT": count,
                "QGRS_SEQUENCE": sequence,
                "QGRS_GSCORE": gscore,
                "QGRS_TETRADS": tetrads,
                "QGRS_LENGTH": length,
            }
        )
    for signature, count in sorted((ref_counter - mut_counter).items()):
        sequence, gscore, tetrads, length = signature
        delta_rows.append(
            {
                "MUTANT_HEADER": mutant_header,
                "REFERENCE_HEADER": reference_header,
                "GENE": manifest_row["GENE"],
                "TRANSCRIPT": manifest_row["TRANSCRIPT"],
                "STRAND": manifest_row["STRAND"],
                "REGION": manifest_row["REGION"],
                "VARIANT_CHROM": manifest_row["VARIANT_CHROM"],
                "VARIANT_POS_1BASED": manifest_row["VARIANT_POS_1BASED"],
                "REF_ALLELE_GENOMIC": manifest_row["REF_ALLELE_GENOMIC"],
                "ALT_ALLELE_GENOMIC": manifest_row["ALT_ALLELE_GENOMIC"],
                "REF_ALLELE_TRANSCRIPT": manifest_row["REF_ALLELE_TRANSCRIPT"],
                "ALT_ALLELE_TRANSCRIPT": manifest_row["ALT_ALLELE_TRANSCRIPT"],
                "TRANSCRIPT_POS1": manifest_row["TRANSCRIPT_POS1"],
                "G_CHANGE_CLASS": manifest_row["G_CHANGE_CLASS"],
                "DELTA_TYPE": "lost",
                "DELTA_COUNT": count,
                "QGRS_SEQUENCE": sequence,
                "QGRS_GSCORE": gscore,
                "QGRS_TETRADS": tetrads,
                "QGRS_LENGTH": length,
            }
        )
    return delta_rows


def classify_effect(ref_count, mut_count, delta_count, delta_max_gscore):
    if ref_count == 0 and mut_count > 0:
        return "qgrs_gained_from_none"
    if ref_count > 0 and mut_count == 0:
        return "qgrs_lost_all"
    if delta_count > 0 and delta_max_gscore > 0:
        return "more_qgrs_higher_max_gscore"
    if delta_count < 0 and delta_max_gscore < 0:
        return "fewer_qgrs_lower_max_gscore"
    if delta_max_gscore > 0:
        return "higher_max_gscore"
    if delta_max_gscore < 0:
        return "lower_max_gscore"
    if delta_count > 0:
        return "more_qgrs_same_max_gscore"
    if delta_count < 0:
        return "fewer_qgrs_same_max_gscore"
    return "no_summary_change"


def write_rows(rows, output_path, fieldnames):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Compare reference QGRS outputs to mutant QGRS outputs.")
    parser.add_argument("--reference-summary", required=True, type=Path)
    parser.add_argument("--reference-hits", required=True, type=Path)
    parser.add_argument("--mutant-summary", required=True, type=Path)
    parser.add_argument("--mutant-hits", required=True, type=Path)
    parser.add_argument("--mutant-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    ref_summary_rows = load_tsv(args.reference_summary)
    ref_hit_rows = load_tsv(args.reference_hits)
    mut_summary_rows = load_tsv(args.mutant_summary)
    mut_hit_rows = load_tsv(args.mutant_hits)
    manifest_rows = load_tsv(args.mutant_manifest)

    ref_summary_by_header = {row["HEADER"]: row for row in ref_summary_rows}
    mut_summary_by_header = {row["HEADER"]: row for row in mut_summary_rows}

    ref_hits_by_header = defaultdict(list)
    for row in ref_hit_rows:
        ref_hits_by_header[row["HEADER"]].append(row)

    mut_hits_by_header = defaultdict(list)
    for row in mut_hit_rows:
        mut_hits_by_header[row["HEADER"]].append(row)

    summary_rows = []
    motif_delta_rows = []

    for manifest_row in manifest_rows:
        reference_header = manifest_row["REFERENCE_HEADER"]
        mutant_header = manifest_row["MUTANT_HEADER"]

        ref_summary = ref_summary_by_header[reference_header]
        mut_summary = mut_summary_by_header[mutant_header]
        ref_hits = ref_hits_by_header[reference_header]
        mut_hits = mut_hits_by_header[mutant_header]

        ref_count = safe_int(ref_summary["QGRS_COUNT"])
        mut_count = safe_int(mut_summary["QGRS_COUNT"])
        ref_max_gscore = safe_int(ref_summary["MAX_GSCORE"])
        mut_max_gscore = safe_int(mut_summary["MAX_GSCORE"])
        ref_max_tetrads = safe_int(ref_summary["MAX_TETRADS"])
        mut_max_tetrads = safe_int(mut_summary["MAX_TETRADS"])

        ref_best = best_hit(ref_hits)
        mut_best = best_hit(mut_hits)
        ref_counter = Counter(motif_signature(row) for row in ref_hits)
        mut_counter = Counter(motif_signature(row) for row in mut_hits)

        gained_count = sum((mut_counter - ref_counter).values())
        lost_count = sum((ref_counter - mut_counter).values())
        shared_count = sum((mut_counter & ref_counter).values())

        delta_count = mut_count - ref_count
        delta_max_gscore = mut_max_gscore - ref_max_gscore
        delta_max_tetrads = mut_max_tetrads - ref_max_tetrads

        summary_rows.append(
            {
                "MUTANT_HEADER": mutant_header,
                "REFERENCE_HEADER": reference_header,
                "GENE": manifest_row["GENE"],
                "GENE_NO": manifest_row["GENE_NO"],
                "TRANSCRIPT": manifest_row["TRANSCRIPT"],
                "STRAND": manifest_row["STRAND"],
                "REGION": manifest_row["REGION"],
                "SEQUENCE_ID": manifest_row["SEQUENCE_ID"],
                "VARIANT_CHROM": manifest_row["VARIANT_CHROM"],
                "VARIANT_POS_1BASED": manifest_row["VARIANT_POS_1BASED"],
                "VARIANT_ID": manifest_row["VARIANT_ID"],
                "ALT_INDEX": manifest_row["ALT_INDEX"],
                "REF_ALLELE_GENOMIC": manifest_row["REF_ALLELE_GENOMIC"],
                "ALT_ALLELE_GENOMIC": manifest_row["ALT_ALLELE_GENOMIC"],
                "REF_ALLELE_TRANSCRIPT": manifest_row["REF_ALLELE_TRANSCRIPT"],
                "ALT_ALLELE_TRANSCRIPT": manifest_row["ALT_ALLELE_TRANSCRIPT"],
                "TRANSCRIPT_POS1": manifest_row["TRANSCRIPT_POS1"],
                "G_CHANGE_CLASS": manifest_row["G_CHANGE_CLASS"],
                "QGRS_OVERLAP_COUNT": manifest_row["QGRS_OVERLAP_COUNT"],
                "QGRS_OVERLAP_IDS": manifest_row["QGRS_OVERLAP_IDS"],
                "REFERENCE_QGRS_COUNT": ref_count,
                "MUTANT_QGRS_COUNT": mut_count,
                "DELTA_QGRS_COUNT": delta_count,
                "REFERENCE_MAX_GSCORE": ref_max_gscore,
                "MUTANT_MAX_GSCORE": mut_max_gscore,
                "DELTA_MAX_GSCORE": delta_max_gscore,
                "REFERENCE_MAX_TETRADS": ref_max_tetrads,
                "MUTANT_MAX_TETRADS": mut_max_tetrads,
                "DELTA_MAX_TETRADS": delta_max_tetrads,
                "REFERENCE_HAS_QGRS": ref_summary["HAS_QGRS"],
                "MUTANT_HAS_QGRS": mut_summary["HAS_QGRS"],
                "GAINED_MOTIF_INSTANCES": gained_count,
                "LOST_MOTIF_INSTANCES": lost_count,
                "SHARED_MOTIF_INSTANCES": shared_count,
                "REFERENCE_TOP_SEQUENCE": ref_best["QGRS_SEQUENCE"] if ref_best else "",
                "MUTANT_TOP_SEQUENCE": mut_best["QGRS_SEQUENCE"] if mut_best else "",
                "REFERENCE_TOP_GSCORE": ref_best["GSCORE"] if ref_best else "",
                "MUTANT_TOP_GSCORE": mut_best["GSCORE"] if mut_best else "",
                "REFERENCE_TOP_TETRADS": ref_best["TETRADS"] if ref_best else "",
                "MUTANT_TOP_TETRADS": mut_best["TETRADS"] if mut_best else "",
                "TOP_SEQUENCE_CHANGED": "y"
                if (ref_best["QGRS_SEQUENCE"] if ref_best else "") != (mut_best["QGRS_SEQUENCE"] if mut_best else "")
                else "n",
                "EFFECT_CLASS": classify_effect(ref_count, mut_count, delta_count, delta_max_gscore),
                "CLNSIG": manifest_row["CLNSIG"],
                "CLNDN": manifest_row["CLNDN"],
                "RS": manifest_row["RS"],
                "ALLELEID": manifest_row["ALLELEID"],
                "CLNVC": manifest_row["CLNVC"],
            }
        )

        motif_delta_rows.extend(
            build_motif_delta_rows(manifest_row, ref_hits, mut_hits, reference_header, mutant_header)
        )

    comparison_summary_fieldnames = [
        "MUTANT_HEADER",
        "REFERENCE_HEADER",
        "GENE",
        "GENE_NO",
        "TRANSCRIPT",
        "STRAND",
        "REGION",
        "SEQUENCE_ID",
        "VARIANT_CHROM",
        "VARIANT_POS_1BASED",
        "VARIANT_ID",
        "ALT_INDEX",
        "REF_ALLELE_GENOMIC",
        "ALT_ALLELE_GENOMIC",
        "REF_ALLELE_TRANSCRIPT",
        "ALT_ALLELE_TRANSCRIPT",
        "TRANSCRIPT_POS1",
        "G_CHANGE_CLASS",
        "QGRS_OVERLAP_COUNT",
        "QGRS_OVERLAP_IDS",
        "REFERENCE_QGRS_COUNT",
        "MUTANT_QGRS_COUNT",
        "DELTA_QGRS_COUNT",
        "REFERENCE_MAX_GSCORE",
        "MUTANT_MAX_GSCORE",
        "DELTA_MAX_GSCORE",
        "REFERENCE_MAX_TETRADS",
        "MUTANT_MAX_TETRADS",
        "DELTA_MAX_TETRADS",
        "REFERENCE_HAS_QGRS",
        "MUTANT_HAS_QGRS",
        "GAINED_MOTIF_INSTANCES",
        "LOST_MOTIF_INSTANCES",
        "SHARED_MOTIF_INSTANCES",
        "REFERENCE_TOP_SEQUENCE",
        "MUTANT_TOP_SEQUENCE",
        "REFERENCE_TOP_GSCORE",
        "MUTANT_TOP_GSCORE",
        "REFERENCE_TOP_TETRADS",
        "MUTANT_TOP_TETRADS",
        "TOP_SEQUENCE_CHANGED",
        "EFFECT_CLASS",
        "CLNSIG",
        "CLNDN",
        "RS",
        "ALLELEID",
        "CLNVC",
    ]
    motif_delta_fieldnames = [
        "MUTANT_HEADER",
        "REFERENCE_HEADER",
        "GENE",
        "TRANSCRIPT",
        "STRAND",
        "REGION",
        "VARIANT_CHROM",
        "VARIANT_POS_1BASED",
        "REF_ALLELE_GENOMIC",
        "ALT_ALLELE_GENOMIC",
        "REF_ALLELE_TRANSCRIPT",
        "ALT_ALLELE_TRANSCRIPT",
        "TRANSCRIPT_POS1",
        "G_CHANGE_CLASS",
        "DELTA_TYPE",
        "DELTA_COUNT",
        "QGRS_SEQUENCE",
        "QGRS_GSCORE",
        "QGRS_TETRADS",
        "QGRS_LENGTH",
    ]

    overall_summary = [
        {
            "TOTAL_MUTANT_COMPARISONS": len(summary_rows),
            "MEAN_DELTA_QGRS_COUNT": round(sum(row["DELTA_QGRS_COUNT"] for row in summary_rows) / len(summary_rows), 3)
            if summary_rows
            else 0.0,
            "MEAN_DELTA_MAX_GSCORE": round(
                sum(row["DELTA_MAX_GSCORE"] for row in summary_rows) / len(summary_rows), 3
            )
            if summary_rows
            else 0.0,
            "SNPS_WITH_HIGHER_MAX_GSCORE": sum(1 for row in summary_rows if row["DELTA_MAX_GSCORE"] > 0),
            "SNPS_WITH_LOWER_MAX_GSCORE": sum(1 for row in summary_rows if row["DELTA_MAX_GSCORE"] < 0),
            "SNPS_WITH_UNCHANGED_MAX_GSCORE": sum(1 for row in summary_rows if row["DELTA_MAX_GSCORE"] == 0),
            "SNPS_WITH_MORE_QGRS": sum(1 for row in summary_rows if row["DELTA_QGRS_COUNT"] > 0),
            "SNPS_WITH_FEWER_QGRS": sum(1 for row in summary_rows if row["DELTA_QGRS_COUNT"] < 0),
            "SNPS_WITH_UNCHANGED_QGRS_COUNT": sum(1 for row in summary_rows if row["DELTA_QGRS_COUNT"] == 0),
            "SNPS_WITH_GAINED_MOTIFS": sum(1 for row in summary_rows if row["GAINED_MOTIF_INSTANCES"] > 0),
            "SNPS_WITH_LOST_MOTIFS": sum(1 for row in summary_rows if row["LOST_MOTIF_INSTANCES"] > 0),
            "SNPS_WITH_TOP_SEQUENCE_CHANGED": sum(1 for row in summary_rows if row["TOP_SEQUENCE_CHANGED"] == "y"),
        }
    ]

    summary_path = args.output_dir / f"{args.prefix}.comparison_summary.tsv"
    motif_delta_path = args.output_dir / f"{args.prefix}.motif_deltas.tsv"
    overall_path = args.output_dir / f"{args.prefix}.overall_summary.tsv"

    write_rows(summary_rows, summary_path, comparison_summary_fieldnames)
    write_rows(motif_delta_rows, motif_delta_path, motif_delta_fieldnames)
    write_rows(overall_summary, overall_path, list(overall_summary[0].keys()))

    print(f"Wrote {summary_path}")
    print(f"Wrote {motif_delta_path}")
    print(f"Wrote {overall_path}")
    print(f"Total mutant comparisons: {len(summary_rows)}")
    print(f"Motif delta rows: {len(motif_delta_rows)}")


if __name__ == "__main__":
    main()
