#!/usr/bin/env python3

# STAR generates SJ.out.tab, which contains nine columns:
# Column    Description
# 1	Chromosome
# 2–3	First and last intronic bases, 1-based
# 4	Strand (0, 1, or 2)
# 5	Splice motif
# 6	Annotated-junction flag
# 7	Uniquely mapping junction-spanning reads
# 8	Multimapping junction-spanning reads
# 9	Maximum overhang

import os
import csv
from collections import defaultdict

CATALOGS = [
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.gencode.catalog.tsv",
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.refseq.catalog.tsv"
]

known = defaultdict(set)
all_rows = []

for path in CATALOGS:
    with open(path) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            all_rows.append(r)
            for pair in r["introns"].split(";"):
                if pair:
                    known[
                        (r["chromosome"], pair, r["strand"])
                    ].add(r["source"])

chrom = all_rows[0]["chromosome"]
start = min(int(r["start"]) for r in all_rows)
end = max(int(r["end"]) for r in all_rows)

samples=[
    "1108_NEU1", "1108_NEU2", "1108_NEU3", "1108_NSCP10",
    "150_N2", "150_NEU_5JAN1", "150_NEU_5JAN2", "150_NSC_P14", "150_NSC_P8",
    "ATC2555_NEU3", "ATC2555_NEU5", "ATC2555_NEU6", "ATC2555_NSCP14", "ATC2555_NSCP8",
    "GOC29_N3", "GOC29_N4"
]

results = []

for sample in samples:

    sj = f"/mnt/faruq2/lab_users/vivek/analysis/RNA-Seq/SCA12/batch2/dragen/{sample}.SJ.out.tab"

    if not os.path.exists(sj):
        raise FileNotFoundError(sj)

    with open(sj) as handle:
        for line in handle:

            v = line.strip().split("\t")

            if len(v) != 9 or v[0] != chrom:
                continue

            s, e = int(v[1]), int(v[2])

            # Keep junctions overlapping the PPP2R2B interval.
            if e < start or s > end:
                continue

            strand_code = int(v[3])
            strand = {0: ".", 1: "+", 2: "-"}[strand_code]

            # Preserve unresolved strand junctions for review.
            if strand == "+":
                continue

            key = (chrom, f"{s}-{e}", strand)

            sources = known.get(key, set())

            results.append({
                "sample": sample,
                "chromosome": chrom,
                "start": s,
                "end": e,
                "strand": strand,
                "motif_code": int(v[4]),
                "known_gencode": int("GENCODE" in sources),
                "known_refseq": int("REFSEQ" in sources),
                "unique_reads": int(v[6]),
                "multi_reads": int(v[7]),
                "max_overhang": int(v[8])
            })

out = "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/results/PPP2R2B_junctions_by_sample.tsv"

with open(out, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=[
            "sample", "chromosome", "start", "end",
            "strand", "motif_code",
            "known_gencode", "known_refseq",
            "unique_reads", "multi_reads", "max_overhang"
        ], delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(results)

print(f"Wrote {len(results)} sample-junction records")
print(out)

#### Summarize novel junction support across samples ####
import pandas as pd

df = pd.read_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/results/PPP2R2B_junctions_by_sample.tsv",
    sep="\t"
)

# Keep sense-strand PPP2R2B junctions.
df = df[df["strand"] == "-"].copy()

summary = (
    df.groupby(
        [
            "chromosome", "start", "end",
            "known_gencode", "known_refseq"
        ],
        as_index=False
    )
    .agg(
        total_unique_reads=("unique_reads", "sum"),
        max_overhang=("max_overhang", "max"),
        sample_count=("sample", "nunique"),
        samples_with_3_reads=(
            "unique_reads",
            lambda x: int((x >= 3).sum())
        )
    )
)

summary["novel_in_both"] = (
    (summary["known_gencode"] == 0)
    & (summary["known_refseq"] == 0)
)

summary.to_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/results/PPP2R2B_junction_summary.tsv",
    sep="\t", index=False
)

# Exploratory priority set, not a definitive validation filter.
candidates = summary[
    summary["novel_in_both"]
    & (summary["total_unique_reads"] >= 5)
    & (summary["max_overhang"] >= 10)
]

candidates.to_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/results/PPP2R2B_prioritized_novel_junctions.tsv",
    sep="\t", index=False
)

print(candidates.to_string(index=False))
