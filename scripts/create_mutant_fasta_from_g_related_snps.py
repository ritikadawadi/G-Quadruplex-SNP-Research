#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_fasta(path):
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


def parse_header(header):
    parsed = {"HEADER": header}
    for item in header.split("|"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def unique_variant_key(row):
    return (
        row["SEQUENCE_ID"],
        row["VARIANT_CHROM"],
        row["VARIANT_POS_1BASED"],
        row["ALT_INDEX"],
        row["REF_ALLELE_GENOMIC"],
        row["ALT_ALLELE_GENOMIC"],
        row["TRANSCRIPT_POS1"],
        row["REF_ALLELE_TRANSCRIPT"],
        row["ALT_ALLELE_TRANSCRIPT"],
    )


def build_mutant_header(base_header, row):
    variant_tag = (
        f"{row['VARIANT_CHROM']}:{row['VARIANT_POS_1BASED']}_"
        f"{row['REF_ALLELE_GENOMIC']}>{row['ALT_ALLELE_GENOMIC']}"
    )
    tx_change = f"{row['REF_ALLELE_TRANSCRIPT']}>{row['ALT_ALLELE_TRANSCRIPT']}"
    extra = [
        f"VARIANT={variant_tag}",
        f"ALT_INDEX={row['ALT_INDEX']}",
        f"TX_POS={row['TRANSCRIPT_POS1']}",
        f"TX_CHANGE={tx_change}",
        f"G_CHANGE_CLASS={row['G_CHANGE_CLASS']}",
    ]
    return base_header + "|" + "|".join(extra)


def build_reference_header_from_row(row):
    return (
        f"GENE={row['GENE']}|GENE_NO={row['GENE_NO']}|TRANSCRIPT={row['TRANSCRIPT']}|"
        f"STRAND={row['STRAND']}|COSMIC={row['COSMIC']}|REGION={row['REGION']}"
    )


def main():
    parser = argparse.ArgumentParser(description="Create one mutant FASTA sequence per unique G-related SNP.")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--g-related-snps", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    args = parser.parse_args()

    reference_records = parse_fasta(args.reference_fasta)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)

    grouped = {}
    qgrs_groups = defaultdict(list)

    with args.g_related_snps.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = unique_variant_key(row)
            grouped.setdefault(key, row)
            qgrs_groups[key].append(row["QGRS_ID"])

    manifest_rows = []
    with args.output_fasta.open("w") as fasta_handle:
        for key in sorted(grouped, key=lambda item: (int(item[0]), item[1], int(item[2]), int(item[3]))):
            row = grouped[key]
            header = build_reference_header_from_row(row)
            ref_seq = reference_records[header]
            tx_pos1 = int(row["TRANSCRIPT_POS1"])
            tx_pos0 = tx_pos1 - 1
            ref_base = row["REF_ALLELE_TRANSCRIPT"]
            alt_base = row["ALT_ALLELE_TRANSCRIPT"]

            observed_base = ref_seq[tx_pos0]
            if observed_base != ref_base:
                raise ValueError(
                    f"Reference base mismatch for {header} at transcript position {tx_pos1}: "
                    f"expected {ref_base}, found {observed_base}"
                )

            mutant_seq = ref_seq[:tx_pos0] + alt_base + ref_seq[tx_pos0 + 1 :]
            mutant_header = build_mutant_header(header, row)
            fasta_handle.write(f">{mutant_header}\n{mutant_seq}\n")

            qgrs_ids = sorted({int(x) for x in qgrs_groups[key]})
            manifest_rows.append(
                {
                    "MUTANT_HEADER": mutant_header,
                    "REFERENCE_HEADER": header,
                    "GENE": row["GENE"],
                    "GENE_NO": row["GENE_NO"],
                    "TRANSCRIPT": row["TRANSCRIPT"],
                    "STRAND": row["STRAND"],
                    "REGION": row["REGION"],
                    "SEQUENCE_ID": row["SEQUENCE_ID"],
                    "VARIANT_CHROM": row["VARIANT_CHROM"],
                    "VARIANT_POS_1BASED": row["VARIANT_POS_1BASED"],
                    "VARIANT_ID": row["VARIANT_ID"],
                    "ALT_INDEX": row["ALT_INDEX"],
                    "REF_ALLELE_GENOMIC": row["REF_ALLELE_GENOMIC"],
                    "ALT_ALLELE_GENOMIC": row["ALT_ALLELE_GENOMIC"],
                    "REF_ALLELE_TRANSCRIPT": ref_base,
                    "ALT_ALLELE_TRANSCRIPT": alt_base,
                    "TRANSCRIPT_POS1": tx_pos1,
                    "REFERENCE_BASE_CONFIRMED": "y",
                    "G_CHANGE_CLASS": row["G_CHANGE_CLASS"],
                    "QGRS_OVERLAP_COUNT": len(qgrs_ids),
                    "QGRS_OVERLAP_IDS": ",".join(str(x) for x in qgrs_ids),
                    "CLNSIG": row["CLNSIG"],
                    "CLNDN": row["CLNDN"],
                    "RS": row["RS"],
                    "ALLELEID": row["ALLELEID"],
                    "CLNVC": row["CLNVC"],
                }
            )

    fieldnames = [
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
        "REFERENCE_BASE_CONFIRMED",
        "G_CHANGE_CLASS",
        "QGRS_OVERLAP_COUNT",
        "QGRS_OVERLAP_IDS",
        "CLNSIG",
        "CLNDN",
        "RS",
        "ALLELEID",
        "CLNVC",
    ]
    with args.output_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Unique mutant sequences written: {len(manifest_rows)}")
    print(f"Wrote {args.output_fasta}")
    print(f"Wrote {args.output_manifest}")


if __name__ == "__main__":
    main()
