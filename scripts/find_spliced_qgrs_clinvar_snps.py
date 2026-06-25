#!/usr/bin/env python3
import argparse
import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


CLINVAR_VCF = Path("/Users/ritika/Desktop/pilot/data/clinvar/clinvar.vcf.gz")


def complement(base):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return base.translate(table).upper()


def chrom_without_chr(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def parse_info(info_text):
    parsed = {}
    for field in info_text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            parsed[key] = value
        else:
            parsed[field] = True
    return parsed


def compress_positions_to_intervals(positions_1based):
    if not positions_1based:
        return []
    sorted_positions = sorted(positions_1based)
    intervals = []
    start = sorted_positions[0]
    end = start
    for pos in sorted_positions[1:]:
        if pos == end + 1:
            end = pos
        else:
            intervals.append((start, end))
            start = pos
            end = pos
    intervals.append((start, end))
    return intervals


def load_bed12_regions(region_bed12):
    regions = {}
    with region_bed12.open() as handle:
        for line in handle:
            chrom, chrom_start, chrom_end, name, score, strand, thick_start, thick_end, rgb, block_count, block_sizes, block_starts = line.rstrip().split("\t")
            chrom_start = int(chrom_start)
            sizes = [int(x) for x in block_sizes.rstrip(",").split(",") if x]
            starts = [int(x) for x in block_starts.rstrip(",").split(",") if x]
            blocks = []
            for size, rel_start in zip(sizes, starts):
                block_start0 = chrom_start + rel_start
                block_end0 = block_start0 + size
                blocks.append((block_start0, block_end0))

            transcript_positions = []
            if strand == "+":
                block_iter = blocks
                for block_start0, block_end0 in block_iter:
                    transcript_positions.extend(pos0 + 1 for pos0 in range(block_start0, block_end0))
            else:
                block_iter = reversed(blocks)
                for block_start0, block_end0 in block_iter:
                    transcript_positions.extend(pos0 + 1 for pos0 in range(block_end0 - 1, block_start0 - 1, -1))

            header = f"{name}({strand})"
            pos_to_txpos1 = {pos1: idx + 1 for idx, pos1 in enumerate(transcript_positions)}
            regions[header] = {
                "chrom": chrom,
                "chrom_query": chrom_without_chr(chrom),
                "name": name,
                "strand": strand,
                "blocks": blocks,
                "tx_pos_to_genomic1": transcript_positions,
                "pos1_to_txpos1": pos_to_txpos1,
                "sequence_length": len(transcript_positions),
            }
    return regions


def load_qgrs_regions(qgrs_hits_path, regions):
    qgrs_rows = []
    by_chrom = defaultdict(list)

    with qgrs_hits_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("STATUS") and row["STATUS"] != "ok":
                continue

            region = regions[row["HEADER"]]
            qgrs_start0 = int(row["START"])
            qgrs_length = int(row["LENGTH"])
            genomic_positions = region["tx_pos_to_genomic1"][qgrs_start0 : qgrs_start0 + qgrs_length]
            genomic_intervals = compress_positions_to_intervals(genomic_positions)
            motif_pos_to_offset = {pos1: idx for idx, pos1 in enumerate(genomic_positions)}

            qgrs_row = {
                "SEQUENCE_ID": row["SEQUENCE_ID"],
                "QGRS_ID": row["QGRS_ID"],
                "HEADER": row["HEADER"],
                "GENE": row["GENE"],
                "GENE_NO": row["GENE_NO"],
                "TRANSCRIPT": row["TRANSCRIPT"],
                "STRAND": row["STRAND"],
                "COSMIC": row["COSMIC"],
                "REGION": row["REGION"],
                "QGRS_START_IN_SEQUENCE": qgrs_start0,
                "QGRS_LENGTH": qgrs_length,
                "QGRS_SEQUENCE": row["QGRS_SEQUENCE"],
                "QGRS_GSCORE": row["GSCORE"],
                "QGRS_TETRADS": row["TETRADS"],
                "CHROM": region["chrom"],
                "CHROM_QUERY": region["chrom_query"],
                "GENOMIC_START_1BASED": min(genomic_positions),
                "GENOMIC_END_1BASED": max(genomic_positions),
                "GENOMIC_INTERVALS": ",".join(f"{start}-{end}" for start, end in genomic_intervals),
                "INTERVAL_COUNT": len(genomic_intervals),
                "motif_pos_to_offset": motif_pos_to_offset,
                "region_pos1_to_txpos1": region["pos1_to_txpos1"],
                "genomic_intervals": genomic_intervals,
            }
            qgrs_rows.append(qgrs_row)
            by_chrom[region["chrom_query"]].append(qgrs_row)

    return qgrs_rows, by_chrom


def write_qgrs_regions(qgrs_rows, output_dir, prefix):
    tsv_path = output_dir / f"{prefix}.qgrs_regions.tsv"
    bed_path = output_dir / f"{prefix}.qgrs_regions.bed"

    fieldnames = [
        "SEQUENCE_ID",
        "QGRS_ID",
        "HEADER",
        "GENE",
        "GENE_NO",
        "TRANSCRIPT",
        "STRAND",
        "COSMIC",
        "REGION",
        "QGRS_START_IN_SEQUENCE",
        "QGRS_LENGTH",
        "QGRS_SEQUENCE",
        "QGRS_GSCORE",
        "QGRS_TETRADS",
        "CHROM",
        "CHROM_QUERY",
        "GENOMIC_START_1BASED",
        "GENOMIC_END_1BASED",
        "GENOMIC_INTERVALS",
        "INTERVAL_COUNT",
    ]

    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in qgrs_rows:
            writer.writerow({key: row[key] for key in fieldnames})

    with bed_path.open("w") as handle:
        for row in qgrs_rows:
            name = f"{row['GENE']}|SEQ={row['SEQUENCE_ID']}|QGRS={row['QGRS_ID']}|GSCORE={row['QGRS_GSCORE']}"
            for index, (start1, end1) in enumerate(row["genomic_intervals"], start=1):
                handle.write(
                    "\t".join(
                        [
                            row["CHROM"],
                            str(start1 - 1),
                            str(end1),
                            f"{name}|BLOCK={index}",
                            "0",
                            row["STRAND"],
                        ]
                    )
                    + "\n"
                )

    return tsv_path, bed_path


def query_clinvar(qgrs_rows):
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as region_handle:
        for row in qgrs_rows:
            for start1, end1 in row["genomic_intervals"]:
                region_handle.write(
                    "\t".join(
                        [
                            row["CHROM_QUERY"],
                            str(start1 - 1),
                            str(end1),
                        ]
                    )
                    + "\n"
                )
        temp_path = Path(region_handle.name)

    try:
        result = subprocess.run(
            ["bcftools", "view", "-R", str(temp_path), "-H", str(CLINVAR_VCF)],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        temp_path.unlink(missing_ok=True)

    return [line for line in result.stdout.splitlines() if line.strip()]


def build_overlap_rows(variant_lines, regions_by_chrom):
    all_rows = []
    g_rows = []

    for line in variant_lines:
        chrom, pos, variant_id, ref, alts, qual, filt, info, *rest = line.rstrip().split("\t")
        pos1 = int(pos)
        info_map = parse_info(info)

        for alt_index, alt in enumerate(alts.split(","), start=1):
            if len(ref) != 1 or len(alt) != 1:
                continue

            for region in regions_by_chrom.get(chrom, []):
                if pos1 not in region["motif_pos_to_offset"]:
                    continue

                qgrs_offset0 = region["motif_pos_to_offset"][pos1]
                qgrs_base = region["QGRS_SEQUENCE"][qgrs_offset0]
                transcript_pos1 = region["region_pos1_to_txpos1"][pos1]

                if region["STRAND"] == "+":
                    ref_tx = ref.upper()
                    alt_tx = alt.upper()
                else:
                    ref_tx = complement(ref)
                    alt_tx = complement(alt)

                if ref_tx == "G" and alt_tx != "G":
                    g_change = "from_G"
                elif ref_tx != "G" and alt_tx == "G":
                    g_change = "to_G"
                elif ref_tx == "G" and alt_tx == "G":
                    g_change = "G_to_G"
                else:
                    g_change = "non_G_change"

                row = {
                    "SEQUENCE_ID": region["SEQUENCE_ID"],
                    "QGRS_ID": region["QGRS_ID"],
                    "GENE": region["GENE"],
                    "GENE_NO": region["GENE_NO"],
                    "TRANSCRIPT": region["TRANSCRIPT"],
                    "STRAND": region["STRAND"],
                    "COSMIC": region["COSMIC"],
                    "REGION": region["REGION"],
                    "CHROM": region["CHROM"],
                    "QGRS_GENOMIC_START_1BASED": region["GENOMIC_START_1BASED"],
                    "QGRS_GENOMIC_END_1BASED": region["GENOMIC_END_1BASED"],
                    "GENOMIC_INTERVALS": region["GENOMIC_INTERVALS"],
                    "QGRS_START_IN_SEQUENCE": region["QGRS_START_IN_SEQUENCE"],
                    "QGRS_LENGTH": region["QGRS_LENGTH"],
                    "QGRS_GSCORE": region["QGRS_GSCORE"],
                    "QGRS_TETRADS": region["QGRS_TETRADS"],
                    "QGRS_SEQUENCE": region["QGRS_SEQUENCE"],
                    "VARIANT_CHROM": chrom,
                    "VARIANT_POS_1BASED": pos1,
                    "VARIANT_ID": variant_id,
                    "ALT_INDEX": alt_index,
                    "REF_ALLELE_GENOMIC": ref.upper(),
                    "ALT_ALLELE_GENOMIC": alt.upper(),
                    "REF_ALLELE_TRANSCRIPT": ref_tx,
                    "ALT_ALLELE_TRANSCRIPT": alt_tx,
                    "TRANSCRIPT_POS1": transcript_pos1,
                    "QGRS_OFFSET1": qgrs_offset0 + 1,
                    "QGRS_BASE_AT_VARIANT": qgrs_base,
                    "REF_MATCHES_QGRS_BASE": "y" if qgrs_base == ref_tx else "n",
                    "G_CHANGE_CLASS": g_change,
                    "IS_G_RELATED": "y" if g_change in {"from_G", "to_G"} else "n",
                    "CLNSIG": info_map.get("CLNSIG", ""),
                    "CLNDN": info_map.get("CLNDN", ""),
                    "RS": info_map.get("RS", ""),
                    "ALLELEID": info_map.get("ALLELEID", ""),
                    "CLNVC": info_map.get("CLNVC", ""),
                    "GENEINFO": info_map.get("GENEINFO", ""),
                }
                all_rows.append(row)
                if row["IS_G_RELATED"] == "y":
                    g_rows.append(row)

    return all_rows, g_rows


def write_overlap_rows(rows, output_path):
    fieldnames = [
        "SEQUENCE_ID",
        "QGRS_ID",
        "GENE",
        "GENE_NO",
        "TRANSCRIPT",
        "STRAND",
        "COSMIC",
        "REGION",
        "CHROM",
        "QGRS_GENOMIC_START_1BASED",
        "QGRS_GENOMIC_END_1BASED",
        "GENOMIC_INTERVALS",
        "QGRS_START_IN_SEQUENCE",
        "QGRS_LENGTH",
        "QGRS_GSCORE",
        "QGRS_TETRADS",
        "QGRS_SEQUENCE",
        "VARIANT_CHROM",
        "VARIANT_POS_1BASED",
        "VARIANT_ID",
        "ALT_INDEX",
        "REF_ALLELE_GENOMIC",
        "ALT_ALLELE_GENOMIC",
        "REF_ALLELE_TRANSCRIPT",
        "ALT_ALLELE_TRANSCRIPT",
        "TRANSCRIPT_POS1",
        "QGRS_OFFSET1",
        "QGRS_BASE_AT_VARIANT",
        "REF_MATCHES_QGRS_BASE",
        "G_CHANGE_CLASS",
        "IS_G_RELATED",
        "CLNSIG",
        "CLNDN",
        "RS",
        "ALLELEID",
        "CLNVC",
        "GENEINFO",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Find ClinVar SNPs in spliced QGRS regions such as 5'UTR or 3'UTR.")
    parser.add_argument("--region-bed12", required=True, type=Path)
    parser.add_argument("--qgrs-hits", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    regions = load_bed12_regions(args.region_bed12)
    qgrs_rows, regions_by_chrom = load_qgrs_regions(args.qgrs_hits, regions)
    qgrs_tsv, qgrs_bed = write_qgrs_regions(qgrs_rows, args.output_dir, args.prefix)

    variant_lines = query_clinvar(qgrs_rows)
    all_rows, g_rows = build_overlap_rows(variant_lines, regions_by_chrom)

    all_path = args.output_dir / f"{args.prefix}.clinvar_snps.tsv"
    g_path = args.output_dir / f"{args.prefix}.clinvar_g_related_snps.tsv"
    write_overlap_rows(all_rows, all_path)
    write_overlap_rows(g_rows, g_path)

    print(f"Wrote {qgrs_tsv}")
    print(f"Wrote {qgrs_bed}")
    print(f"Wrote {all_path}")
    print(f"Wrote {g_path}")
    print(f"Total ClinVar SNP overlaps: {len(all_rows)}")
    print(f"G-related overlaps: {len(g_rows)}")


if __name__ == "__main__":
    main()
