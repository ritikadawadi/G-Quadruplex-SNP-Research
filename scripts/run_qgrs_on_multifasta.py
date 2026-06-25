#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


QGRS_BIN = Path("/Users/ritika/Desktop/pilot/qgrs-cpp/bin/qgrs")


def parse_fasta(path):
    header = None
    seq_parts = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def parse_header(header):
    parsed = {"HEADER": header}
    for item in header.split("|"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def run_qgrs(sequence, timeout_seconds):
    with tempfile.NamedTemporaryFile("w", delete=False) as temp_handle:
        temp_handle.write(sequence + "\n")
        temp_path = temp_handle.name

    command = [str(QGRS_BIN), "-json", "-i", temp_path]
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        Path(temp_path).unlink(missing_ok=True)
        return None, "timeout"
    Path(temp_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"qgrs failed with exit code {result.returncode}: {result.stderr.strip()}")

    stdout = result.stdout.strip()
    if not stdout:
        return [], "ok"
    return json.loads(stdout).get("results", []), "ok"


def process_record(task):
    sequence_id, header, sequence, timeout_seconds = task
    header_info = parse_header(header)
    qgrs_hits, status = run_qgrs(sequence, timeout_seconds)

    if status == "timeout":
        summary_row = {
            "SEQUENCE_ID": sequence_id,
            "HEADER": header,
            "GENE": header_info.get("GENE", ""),
            "GENE_NO": header_info.get("GENE_NO", ""),
            "TRANSCRIPT": header_info.get("TRANSCRIPT", ""),
            "STRAND": header_info.get("STRAND", ""),
            "COSMIC": header_info.get("COSMIC", ""),
            "REGION": header_info.get("REGION", ""),
            "INPUT_SEQUENCE_LENGTH": len(sequence),
            "QGRS_COUNT": "",
            "MAX_GSCORE": "",
            "MAX_TETRADS": "",
            "HAS_QGRS": "timeout",
            "STATUS": "timeout",
        }
        return [], summary_row

    hit_rows = []
    for qgrs_id, hit in enumerate(qgrs_hits, start=1):
        hit_rows.append(
            {
                "SEQUENCE_ID": sequence_id,
                "HEADER": header,
                "GENE": header_info.get("GENE", ""),
                "GENE_NO": header_info.get("GENE_NO", ""),
                "TRANSCRIPT": header_info.get("TRANSCRIPT", ""),
                "STRAND": header_info.get("STRAND", ""),
                "COSMIC": header_info.get("COSMIC", ""),
                "REGION": header_info.get("REGION", ""),
                "QGRS_ID": qgrs_id,
                "START": hit.get("start", ""),
                "TETRAD1": hit.get("tetrad1", ""),
                "TETRAD2": hit.get("tetrad2", ""),
                "TETRAD3": hit.get("tetrad3", ""),
                "TETRAD4": hit.get("tetrad4", ""),
                "Y1": hit.get("y1", ""),
                "Y2": hit.get("y2", ""),
                "Y3": hit.get("y3", ""),
                "TETRADS": hit.get("tetrads", ""),
                "LENGTH": hit.get("length", ""),
                "GSCORE": hit.get("gscore", ""),
                "QGRS_SEQUENCE": hit.get("sequence", ""),
                "INPUT_SEQUENCE_LENGTH": len(sequence),
                    "STATUS": "ok",
            }
        )

    max_gscore = max((hit.get("gscore", 0) for hit in qgrs_hits), default=0)
    max_tetrads = max((hit.get("tetrads", 0) for hit in qgrs_hits), default=0)
    summary_row = {
        "SEQUENCE_ID": sequence_id,
        "HEADER": header,
        "GENE": header_info.get("GENE", ""),
        "GENE_NO": header_info.get("GENE_NO", ""),
        "TRANSCRIPT": header_info.get("TRANSCRIPT", ""),
        "STRAND": header_info.get("STRAND", ""),
        "COSMIC": header_info.get("COSMIC", ""),
        "REGION": header_info.get("REGION", ""),
        "INPUT_SEQUENCE_LENGTH": len(sequence),
        "QGRS_COUNT": len(qgrs_hits),
        "MAX_GSCORE": max_gscore,
        "MAX_TETRADS": max_tetrads,
        "HAS_QGRS": "y" if qgrs_hits else "n",
        "STATUS": "ok",
    }
    return hit_rows, summary_row


def write_hits(rows, output_path):
    fieldnames = [
        "SEQUENCE_ID",
        "HEADER",
        "GENE",
        "GENE_NO",
        "TRANSCRIPT",
        "STRAND",
        "COSMIC",
        "REGION",
        "QGRS_ID",
        "START",
        "TETRAD1",
        "TETRAD2",
        "TETRAD3",
        "TETRAD4",
        "Y1",
        "Y2",
        "Y3",
        "TETRADS",
        "LENGTH",
        "GSCORE",
        "QGRS_SEQUENCE",
        "INPUT_SEQUENCE_LENGTH",
        "STATUS",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows, output_path):
    fieldnames = [
        "SEQUENCE_ID",
        "HEADER",
        "GENE",
        "GENE_NO",
        "TRANSCRIPT",
        "STRAND",
        "COSMIC",
        "REGION",
        "INPUT_SEQUENCE_LENGTH",
        "QGRS_COUNT",
        "MAX_GSCORE",
        "MAX_TETRADS",
        "HAS_QGRS",
        "STATUS",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Run qgrs on each record in a multi-FASTA file.")
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--per-sequence-timeout", type=int, default=120)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    hit_rows = []
    summary_rows = []
    tasks = [
        (sequence_id, header, sequence, args.per_sequence_timeout)
        for sequence_id, (header, sequence) in enumerate(parse_fasta(args.input_fasta), start=1)
    ]

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        future_to_sequence = {executor.submit(process_record, task): task[0] for task in tasks}
        completed = 0
        for future in as_completed(future_to_sequence):
            record_hit_rows, summary_row = future.result()
            hit_rows.extend(record_hit_rows)
            summary_rows.append(summary_row)
            completed += 1
            if completed % 25 == 0 or completed == len(tasks):
                print(f"Completed {completed}/{len(tasks)} sequences", flush=True)

    hit_rows.sort(key=lambda row: (int(row["SEQUENCE_ID"]), int(row["QGRS_ID"])))
    summary_rows.sort(key=lambda row: int(row["SEQUENCE_ID"]))

    hits_path = args.output_dir / f"{args.prefix}.qgrs_hits.tsv"
    summary_path = args.output_dir / f"{args.prefix}.qgrs_summary.tsv"
    write_hits(hit_rows, hits_path)
    write_summary(summary_rows, summary_path)

    print(f"Processed {len(summary_rows)} FASTA records")
    print(f"Wrote {hits_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
