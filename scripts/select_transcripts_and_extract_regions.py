#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path("/Users/ritika/Desktop/pilot")
GENE_LIST = PROJECT_ROOT / "gene" / "oncogens.tsv"
GTF_PATH = PROJECT_ROOT / "data" / "reference" / "gencode.v46.annotation.gtf"
FAI_PATH = PROJECT_ROOT / "data" / "reference" / "GRCh38.primary_assembly.genome.fa.fai"
REGION_DIR = PROJECT_ROOT / "data" / "regions"


def parse_attributes(attr_text):
    attrs = defaultdict(list)
    for raw_field in attr_text.strip().split(";"):
        field = raw_field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        key = parts[0]
        value = parts[1].strip().strip('"') if len(parts) > 1 else ""
        attrs[key].append(value)
    parsed = {}
    for key, values in attrs.items():
        if key == "tag":
            parsed[key] = values
        else:
            parsed[key] = values[-1]
    return parsed


def strip_version(identifier):
    return identifier.split(".")[0] if identifier else ""


def parse_tsl(value):
    if not value or value == "NA":
        return 99
    try:
        return int(value)
    except ValueError:
        return 99


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_length(intervals):
    return sum(end - start + 1 for start, end in intervals)


def make_bed12_line(chrom, name, strand, intervals):
    merged = merge_intervals(intervals)
    chrom_start = merged[0][0] - 1
    chrom_end = merged[-1][1]
    block_sizes = ",".join(str(end - start + 1) for start, end in merged)
    block_starts = ",".join(str((start - 1) - chrom_start) for start, _ in merged)
    return "\t".join(
        [
            chrom,
            str(chrom_start),
            str(chrom_end),
            name,
            "0",
            strand,
            str(chrom_start),
            str(chrom_end),
            "0",
            str(len(merged)),
            block_sizes,
            block_starts,
        ]
    )


def make_bed6_line(chrom, start0, end0, name, strand):
    return "\t".join([chrom, str(start0), str(end0), name, "0", strand])


def transcript_rank_key(tx):
    tags = set(tx["tags"])
    return (
        1 if tx["transcript_type"] == "protein_coding" else 0,
        1 if tx["five_utr_length"] > 0 and tx["three_utr_length"] > 0 else 0,
        1 if tx["five_utr_length"] > 0 or tx["three_utr_length"] > 0 else 0,
        1 if "basic" in tags else 0,
        1 if "GENCODE_Primary" in tags else 0,
        -parse_tsl(tx["transcript_support_level"]),
        tx["exon_length"],
        tx["transcript_span"],
        tx["transcript_id"],
    )


def choose_transcript(transcripts):
    def pick(candidates, tier_name):
        if not candidates:
            return None, None
        return sorted(candidates, key=transcript_rank_key, reverse=True)[0], tier_name

    mane_with_utr = [tx for tx in transcripts if "MANE_Select" in tx["tags"] and (tx["five_utr_length"] > 0 or tx["three_utr_length"] > 0)]
    selected, tier = pick(mane_with_utr, "MANE_Select_with_UTR")
    if selected:
        return selected, tier

    canonical_with_utr = [tx for tx in transcripts if "Ensembl_canonical" in tx["tags"] and (tx["five_utr_length"] > 0 or tx["three_utr_length"] > 0)]
    selected, tier = pick(canonical_with_utr, "Ensembl_canonical_with_UTR")
    if selected:
        return selected, tier

    protein_with_utr = [tx for tx in transcripts if tx["transcript_type"] == "protein_coding" and (tx["five_utr_length"] > 0 or tx["three_utr_length"] > 0)]
    selected, tier = pick(protein_with_utr, "longest_protein_coding_with_UTR")
    if selected:
        return selected, tier

    mane = [tx for tx in transcripts if "MANE_Select" in tx["tags"]]
    selected, tier = pick(mane, "MANE_Select")
    if selected:
        return selected, tier

    canonical = [tx for tx in transcripts if "Ensembl_canonical" in tx["tags"]]
    selected, tier = pick(canonical, "Ensembl_canonical")
    if selected:
        return selected, tier

    return sorted(transcripts, key=transcript_rank_key, reverse=True)[0], "longest_annotated_transcript"


def load_genes():
    with GENE_LIST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    targets = []
    gene_ids = {}
    gene_symbols = {}
    for row in rows:
        gene_id_base = strip_version((row.get("ENSEMBL_GENE_ID_BASE") or row.get("ENSEMBL_GENE_ID") or "").strip())
        symbol = (row.get("GENE_SYMBOL") or "").strip()
        record = {
            "gene_no": row["GENE_NO"],
            "gene_symbol": symbol,
            "cosmic_gene_id": row["COSMIC_GENE_ID"],
            "ensembl_gene_id_base": gene_id_base,
            "source_row": row,
        }
        targets.append(record)
        if gene_id_base:
            gene_ids[gene_id_base] = record
        if symbol:
            gene_symbols[symbol] = record
    return targets, gene_ids, gene_symbols


def load_chrom_sizes():
    sizes = {}
    with FAI_PATH.open() as handle:
        for line in handle:
            chrom, length, *_ = line.rstrip().split("\t")
            sizes[chrom] = int(length)
    return sizes


def build_transcript_store(target_gene_ids, target_symbols):
    transcript_store = {}
    with GTF_PATH.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_text = line.rstrip().split("\t")
            attrs = parse_attributes(attrs_text)
            gene_id_base = strip_version(attrs.get("gene_id", ""))
            gene_name = attrs.get("gene_name", "")

            matched = None
            match_basis = None
            if gene_id_base and gene_id_base in target_gene_ids:
                matched = target_gene_ids[gene_id_base]
                match_basis = "ensembl_gene_id"
            elif gene_name and gene_name in target_symbols:
                matched = target_symbols[gene_name]
                match_basis = "gene_symbol"
            if matched is None:
                continue

            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            tx = transcript_store.setdefault(
                transcript_id,
                {
                    "gene_no": matched["gene_no"],
                    "gene_symbol": matched["gene_symbol"],
                    "cosmic_gene_id": matched["cosmic_gene_id"],
                    "requested_gene_id": matched["ensembl_gene_id_base"],
                    "matched_by": match_basis,
                    "chrom": chrom,
                    "source": source,
                    "strand": strand,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_id_base": gene_id_base,
                    "gene_name": gene_name,
                    "gene_type": attrs.get("gene_type", ""),
                    "transcript_id": transcript_id,
                    "transcript_name": attrs.get("transcript_name", ""),
                    "transcript_type": attrs.get("transcript_type", ""),
                    "transcript_support_level": attrs.get("transcript_support_level", ""),
                    "tags": list(attrs.get("tag", [])),
                    "start": None,
                    "end": None,
                    "exons": [],
                    "utrs": [],
                    "cds": [],
                },
            )

            start_i = int(start)
            end_i = int(end)
            if feature == "transcript":
                tx["start"] = start_i
                tx["end"] = end_i
                tx["chrom"] = chrom
                tx["strand"] = strand
                tx["gene_type"] = attrs.get("gene_type", tx["gene_type"])
                tx["transcript_name"] = attrs.get("transcript_name", tx["transcript_name"])
                tx["transcript_type"] = attrs.get("transcript_type", tx["transcript_type"])
                tx["transcript_support_level"] = attrs.get("transcript_support_level", tx["transcript_support_level"])
                existing_tags = set(tx["tags"])
                existing_tags.update(attrs.get("tag", []))
                tx["tags"] = sorted(existing_tags)
            elif feature == "exon":
                tx["exons"].append((start_i, end_i))
            elif feature == "UTR":
                tx["utrs"].append((start_i, end_i))
            elif feature == "CDS":
                tx["cds"].append((start_i, end_i))
    return transcript_store


def annotate_transcripts(transcript_store):
    genes = defaultdict(list)
    for tx in transcript_store.values():
        tx["exons"] = merge_intervals(tx["exons"])
        tx["utrs"] = merge_intervals(tx["utrs"])
        tx["cds"] = merge_intervals(tx["cds"])
        tx["exon_length"] = interval_length(tx["exons"])
        tx["transcript_span"] = (tx["end"] - tx["start"] + 1) if tx["start"] and tx["end"] else 0

        five_utrs = []
        three_utrs = []
        if tx["cds"]:
            cds_start = min(start for start, _ in tx["cds"])
            cds_end = max(end for _, end in tx["cds"])
            if tx["strand"] == "+":
                five_utrs = [(start, end) for start, end in tx["utrs"] if end < cds_start]
                three_utrs = [(start, end) for start, end in tx["utrs"] if start > cds_end]
            else:
                five_utrs = [(start, end) for start, end in tx["utrs"] if start > cds_end]
                three_utrs = [(start, end) for start, end in tx["utrs"] if end < cds_start]

        tx["five_utrs"] = merge_intervals(five_utrs)
        tx["three_utrs"] = merge_intervals(three_utrs)
        tx["five_utr_length"] = interval_length(tx["five_utrs"])
        tx["three_utr_length"] = interval_length(tx["three_utrs"])
        genes[tx["gene_no"]].append(tx)
    return genes


def write_selected_transcripts(targets, transcripts_by_gene, chrom_sizes):
    selected_rows = []
    five_utr_lines = []
    three_utr_lines = []
    tss_lines = []

    fieldnames = [
        "GENE_NO",
        "GENE_SYMBOL",
        "COSMIC_GENE_ID",
        "REQUESTED_ENSEMBL_GENE_ID_BASE",
        "MATCHED_GENE_ID",
        "MATCHED_BY",
        "TRANSCRIPT_ID",
        "TRANSCRIPT_NAME",
        "TRANSCRIPT_TYPE",
        "GENE_TYPE",
        "CHROMOSOME",
        "STRAND",
        "TRANSCRIPT_START",
        "TRANSCRIPT_END",
        "TRANSCRIPT_SPAN",
        "EXON_LENGTH",
        "FIVE_UTR_LENGTH",
        "THREE_UTR_LENGTH",
        "HAS_FIVE_UTR",
        "HAS_THREE_UTR",
        "TSS_UPSTREAM500_LENGTH",
        "SELECTION_TIER",
        "TRANSCRIPT_SUPPORT_LEVEL",
        "TAGS",
    ]

    for target in targets:
        gene_no = target["gene_no"]
        transcripts = transcripts_by_gene.get(gene_no, [])
        if not transcripts:
            selected_rows.append(
                {
                    "GENE_NO": gene_no,
                    "GENE_SYMBOL": target["gene_symbol"],
                    "COSMIC_GENE_ID": target["cosmic_gene_id"],
                    "REQUESTED_ENSEMBL_GENE_ID_BASE": target["ensembl_gene_id_base"],
                    "MATCHED_GENE_ID": "",
                    "MATCHED_BY": "",
                    "TRANSCRIPT_ID": "",
                    "TRANSCRIPT_NAME": "",
                    "TRANSCRIPT_TYPE": "",
                    "GENE_TYPE": "",
                    "CHROMOSOME": "",
                    "STRAND": "",
                    "TRANSCRIPT_START": "",
                    "TRANSCRIPT_END": "",
                    "TRANSCRIPT_SPAN": "",
                    "EXON_LENGTH": "",
                    "FIVE_UTR_LENGTH": "",
                    "THREE_UTR_LENGTH": "",
                    "HAS_FIVE_UTR": "n",
                    "HAS_THREE_UTR": "n",
                    "TSS_UPSTREAM500_LENGTH": "",
                    "SELECTION_TIER": "no_transcript_match",
                    "TRANSCRIPT_SUPPORT_LEVEL": "",
                    "TAGS": "",
                }
            )
            continue

        selected, tier = choose_transcript(transcripts)
        chrom_size = chrom_sizes[selected["chrom"]]

        if selected["strand"] == "+":
            tss_start0 = max(0, selected["start"] - 1 - 500)
            tss_end0 = selected["start"] - 1
        else:
            tss_start0 = selected["end"]
            tss_end0 = min(chrom_size, selected["end"] + 500)
        tss_length = max(0, tss_end0 - tss_start0)

        tags_text = ",".join(selected["tags"])
        selected_rows.append(
            {
                "GENE_NO": selected["gene_no"],
                "GENE_SYMBOL": selected["gene_symbol"],
                "COSMIC_GENE_ID": selected["cosmic_gene_id"],
                "REQUESTED_ENSEMBL_GENE_ID_BASE": selected["requested_gene_id"],
                "MATCHED_GENE_ID": selected["gene_id_base"],
                "MATCHED_BY": selected["matched_by"],
                "TRANSCRIPT_ID": selected["transcript_id"],
                "TRANSCRIPT_NAME": selected["transcript_name"],
                "TRANSCRIPT_TYPE": selected["transcript_type"],
                "GENE_TYPE": selected["gene_type"],
                "CHROMOSOME": selected["chrom"],
                "STRAND": selected["strand"],
                "TRANSCRIPT_START": selected["start"],
                "TRANSCRIPT_END": selected["end"],
                "TRANSCRIPT_SPAN": selected["transcript_span"],
                "EXON_LENGTH": selected["exon_length"],
                "FIVE_UTR_LENGTH": selected["five_utr_length"],
                "THREE_UTR_LENGTH": selected["three_utr_length"],
                "HAS_FIVE_UTR": "y" if selected["five_utr_length"] > 0 else "n",
                "HAS_THREE_UTR": "y" if selected["three_utr_length"] > 0 else "n",
                "TSS_UPSTREAM500_LENGTH": tss_length,
                "SELECTION_TIER": tier,
                "TRANSCRIPT_SUPPORT_LEVEL": selected["transcript_support_level"],
                "TAGS": tags_text,
            }
        )

        common_name = (
            f"GENE={selected['gene_symbol']}|GENE_NO={selected['gene_no']}|"
            f"TRANSCRIPT={selected['transcript_id']}|STRAND={selected['strand']}|"
            f"COSMIC={selected['cosmic_gene_id']}"
        )

        if selected["five_utrs"]:
            five_utr_lines.append(
                make_bed12_line(
                    selected["chrom"],
                    common_name + "|REGION=5UTR",
                    selected["strand"],
                    selected["five_utrs"],
                )
            )
        if selected["three_utrs"]:
            three_utr_lines.append(
                make_bed12_line(
                    selected["chrom"],
                    common_name + "|REGION=3UTR",
                    selected["strand"],
                    selected["three_utrs"],
                )
            )
        if tss_length > 0:
            tss_lines.append(
                make_bed6_line(
                    selected["chrom"],
                    tss_start0,
                    tss_end0,
                    common_name + "|REGION=TSS_UPSTREAM500",
                    selected["strand"],
                )
            )

    selected_path = REGION_DIR / "selected_transcripts.tsv"
    with selected_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_rows)

    (REGION_DIR / "5utr.bed12").write_text("\n".join(five_utr_lines) + ("\n" if five_utr_lines else ""))
    (REGION_DIR / "3utr.bed12").write_text("\n".join(three_utr_lines) + ("\n" if three_utr_lines else ""))
    (REGION_DIR / "tss_upstream500.bed").write_text("\n".join(tss_lines) + ("\n" if tss_lines else ""))

    return selected_rows


def main():
    targets, target_gene_ids, target_symbols = load_genes()
    chrom_sizes = load_chrom_sizes()
    transcript_store = build_transcript_store(target_gene_ids, target_symbols)
    transcripts_by_gene = annotate_transcripts(transcript_store)
    selected_rows = write_selected_transcripts(targets, transcripts_by_gene, chrom_sizes)

    matched = sum(1 for row in selected_rows if row["TRANSCRIPT_ID"])
    five_utr = sum(1 for row in selected_rows if row["HAS_FIVE_UTR"] == "y")
    three_utr = sum(1 for row in selected_rows if row["HAS_THREE_UTR"] == "y")
    print(f"Selected transcripts for {matched}/{len(selected_rows)} genes")
    print(f"Genes with 5'UTR: {five_utr}")
    print(f"Genes with 3'UTR: {three_utr}")
    print(f"Wrote {REGION_DIR / 'selected_transcripts.tsv'}")
    print(f"Wrote {REGION_DIR / '5utr.bed12'}")
    print(f"Wrote {REGION_DIR / '3utr.bed12'}")
    print(f"Wrote {REGION_DIR / 'tss_upstream500.bed'}")


if __name__ == "__main__":
    main()
