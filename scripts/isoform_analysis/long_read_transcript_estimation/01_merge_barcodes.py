#!/usr/bin/env python3
import csv
import gzip
from pathlib import Path
from collections import defaultdict

MANIFEST = Path("/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/5RACE_info/samplesheet_5RACE_PCR.csv")
ROOT = Path("/mnt/faruq2/lab_data/raw_data/Nanopore/MF_5RACE_PPP2R2B_29082025/MF_5RACE_PPP2R2B_29082025/20250829_1641_X3_FAY49340_b6e09021/fastq_pass")
OUT = Path("/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read/fastqs")
OUT.mkdir(parents=True, exist_ok=True)

assignments = defaultdict(list)

with MANIFEST.open() as fh:
    for row in csv.DictReader(fh, delimiter=","):
        assignments[row["alias"]].append(row["barcode"])

for sample, barcode_assignments in assignments.items():
    out_fastq = OUT / f"{sample}.fastq.gz"
    total = 0

    with gzip.open(out_fastq, "wt") as dest:
        for barcode in barcode_assignments:
            folder = ROOT / barcode
            files = sorted(folder.glob("*.fastq.gz"))

            if not files:
                raise RuntimeError(
                    f"No FASTQs for {barcode}"
                )

            barcode_reads = 0

            for file in files:
                with gzip.open(file, "rt") as src:
                    while True:
                        header = src.readline()
                        if not header:
                            break

                        sequence = src.readline()
                        plus = src.readline()
                        quality = src.readline()

                        if not sequence or not plus or not quality:
                            raise ValueError(
                                f"Incomplete FASTQ record: {file}"
                            )

                        # Prefix read IDs with the barcode to avoid collisions.
                        original = header[1:].rstrip("\n")
                        new_header = (f"@{barcode}.{original}\n")

                        dest.write(new_header)
                        dest.write(sequence)
                        dest.write(plus)
                        dest.write(quality)

                        barcode_reads += 1

            total += barcode_reads
            print(barcode, sample, barcode_reads, sep="\t")

    print(f"Finished {sample}: {total} reads")
