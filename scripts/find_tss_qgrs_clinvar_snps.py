#!/usr/bin/env python3
import argparse
import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path("/Users/ritika/Desktop/pilot")
TSS_BED = PROJECT_ROOT / "data" / "regions" / "tss_upstream500.bed"
CLINVAR_VCF = PROJECT_ROOT / "data" / "clinvar" / "clinvar.vcf.gz"


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


def load_tss_regions():
    regions = {}
    with TSS_BED.open() as handle:
        for line in handle:
            chrom, start0, end0, name, score, strand = line.rstrip().split("\t")
            header = f"{name}({strand})"
            regions[header] = {
                "chrom": chrom,
                "chrom_query": chrom_without_chr(chrom),
                "start0": int(start0),
                "end0": int(end0),
                "name": name,
                "strand": strand,
            }
    return regions


def load_qgrs_regions(tss_regions):
    region_rows = []
    by_chrom = defaultdict(list)
    with QGRS_HITS.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tss = tss_regions[row["HEADER"]]
            qgrs_start = int(row["START"])
            qgrs_length = int(row["LENGTH"])

            if tss["strand"] == "+":
                genomic_start0 = tss["start0"] + qgrs_start
                genomic_end0 = genomic_start0 + qgrs_length
            else:
                genomic_start0 = tss["end0"] - (qgrs_start + qgrs_length)
                genomic_end0 = tss["end0"] - qgrs_start

            region_row = {
                "SEQUENCE_ID": row["SEQUENCE_ID"],
                "QGRS_ID": row["QGRS_ID"],
                "HEADER": row["HEADER"],
                "GENE": row["GENE"],
                "GENE_NO": row["GENE_NO"],
                "TRANSCRIPT": row["TRANSCRIPT"],
                "STRAND": row["STRAND"],
                "COSMIC": row["COSMIC"],
                "REGION": row["REGION"],
                "QGRS_START_IN_SEQUENCE": qgrs_start,
                "QGRS_LENGTH": qgrs_length,
                "QGRS_SEQUENCE": row["QGRS_SEQUENCE"],
                "QGRS_GSCORE": row["GSCORE"],
                "QGRS_TETRADS": row["TETRADS"],
                "CHROM": tss["chrom"],
                "CHROM_QUERY": tss["chrom_query"],
                "GENOMIC_START0": genomic_start0,
                "GENOMIC_END0": genomic_end0,
                "GENOMIC_START_1BASED": genomic_start0 + 1,
                "GENOMIC_END_1BASED": genomic_end0,
                "TSS_WINDOW_START0": tss["start0"],
                "TSS_WINDOW_END0": tss["end0"],
            }
            region_rows.append(region_row)
            by_chrom[tss["chrom_query"]].append(region_row)

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: (item["GENOMIC_START0"], item["GENOMIC_END0"]))
    return region_rows, by_chrom


def write_qgrs_region_outputs(region_rows, output_dir):
    tsv_path = output_dir / "tss_upstream500_reference.qgrs_regions.tsv"
    bed_path = output_dir / "tss_upstream500_reference.qgrs_regions.bed"

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
        "GENOMIC_START0",
        "GENOMIC_END0",
        "GENOMIC_START_1BASED",
        "GENOMIC_END_1BASED",
        "TSS_WINDOW_START0",
        "TSS_WINDOW_END0",
    ]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(region_rows)

    with bed_path.open("w") as handle:
        for row in region_rows:
            name = f"{row['GENE']}|SEQ={row['SEQUENCE_ID']}|QGRS={row['QGRS_ID']}|GSCORE={row['QGRS_GSCORE']}"
            handle.write(
                "\t".join(
                    [
                        row["CHROM"],
                        str(row["GENOMIC_START0"]),
                        str(row["GENOMIC_END0"]),
                        name,
                        "0",
                        row["STRAND"],
                    ]
                )
                + "\n"
            )

    return tsv_path, bed_path


def query_clinvar(region_rows):
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as region_handle:
        for row in region_rows:
            region_handle.write(
                "\t".join(
                    [
                        row["CHROM_QUERY"],
                        str(row["GENOMIC_START0"]),
                        str(row["GENOMIC_END0"]),
                    ]
                )
                + "\n"
            )
        region_path = Path(region_handle.name)

    try:
        result = subprocess.run(
            ["bcftools", "view", "-R", str(region_path), "-H", str(CLINVAR_VCF)],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        region_path.unlink(missing_ok=True)

    return [line for line in result.stdout.splitlines() if line.strip()]


def build_overlap_rows(variant_lines, regions_by_chrom):
    all_snp_rows = []
    g_related_rows = []

    for line in variant_lines:
        columns = line.rstrip().split("\t")
        chrom, pos, variant_id, ref, alts, qual, filt, info = columns[:8]
        pos1 = int(pos)
        pos0 = pos1 - 1
        info_map = parse_info(info)
        alt_list = alts.split(",")

        for alt_index, alt in enumerate(alt_list, start=1):
            if len(ref) != 1 or len(alt) != 1:
                continue

            for region in regions_by_chrom.get(chrom, []):
                if not (region["GENOMIC_START0"] <= pos0 < region["GENOMIC_END0"]):
                    continue

                if region["STRAND"] == "+":
                    transcript_pos0 = pos0 - region["TSS_WINDOW_START0"]
                else:
                    transcript_pos0 = region["TSS_WINDOW_END0"] - pos1

                qgrs_offset = transcript_pos0 - region["QGRS_START_IN_SEQUENCE"]
                qgrs_base = region["QGRS_SEQUENCE"][qgrs_offset]

                ref_tx = ref.upper() if region["STRAND"] == "+" else complement(ref)
                alt_tx = alt.upper() if region["STRAND"] == "+" else complement(alt)

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
                    "TRANSCRIPT_POS0": transcript_pos0,
                    "TRANSCRIPT_POS1": transcript_pos0 + 1,
                    "QGRS_OFFSET0": qgrs_offset,
                    "QGRS_OFFSET1": qgrs_offset + 1,
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
                all_snp_rows.append(row)
                if row["IS_G_RELATED"] == "y":
                    g_related_rows.append(row)

    return all_snp_rows, g_related_rows


def write_rows(rows, output_path):
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
        "TRANSCRIPT_POS0",
        "TRANSCRIPT_POS1",
        "QGRS_OFFSET0",
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
    parser = argparse.ArgumentParser(description="Find ClinVar SNPs overlapping TSS QGRS regions.")
    parser.add_argument(
        "--qgrs-hits",
        type=Path,
        default=PROJECT_ROOT / "results_tss" / "tss_upstream500_reference.qgrs_hits.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results_tss",
    )
    args = parser.parse_args()

    global QGRS_HITS
    QGRS_HITS = args.qgrs_hits
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tss_regions = load_tss_regions()
    qgrs_region_rows, regions_by_chrom = load_qgrs_regions(tss_regions)
    qgrs_tsv, qgrs_bed = write_qgrs_region_outputs(qgrs_region_rows, args.output_dir)

    variant_lines = query_clinvar(qgrs_region_rows)
    all_snp_rows, g_related_rows = build_overlap_rows(variant_lines, regions_by_chrom)

    all_snps_path = args.output_dir / "tss_upstream500_reference.clinvar_snps.tsv"
    g_related_path = args.output_dir / "tss_upstream500_reference.clinvar_g_related_snps.tsv"
    write_rows(all_snp_rows, all_snps_path)
    write_rows(g_related_rows, g_related_path)

    print(f"Wrote {qgrs_tsv}")
    print(f"Wrote {qgrs_bed}")
    print(f"Wrote {all_snps_path}")
    print(f"Wrote {g_related_path}")
    print(f"Total ClinVar SNP overlaps: {len(all_snp_rows)}")
    print(f"G-related overlaps: {len(g_related_rows)}")


if __name__ == "__main__":
    main()
