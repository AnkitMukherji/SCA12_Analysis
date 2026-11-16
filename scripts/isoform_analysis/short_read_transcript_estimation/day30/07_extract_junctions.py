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
    "102A_S720_L007", "102B_S721_L007", "104A_S761_L008", "104B_S762_L008",
    "105A_S763_L008", "105C_S764_L008", "107A_S765_L008", "107B_S766_L008",
    "109A_S767_L008", "109B_S768_L008", "112A_S769_L008", "112B_S770_L008",
    "113A_S771_L008", "113B_S772_L008", "114A_S773_L008", "114B_S774_L008",
    "115A_S775_L008", "115B_S776_L008", "116A_S777_L008", "116B_S778_L008",
    "117A_S779_L008", "117B_S780_L008", "24A_S713_L007", "24B_S714_L007",
    "25A_S590_L006", "25B_S715_L007", "28A_S716_L007", "28B_S717_L007",
    "29A_S718_L007", "29B_S719_L007"
]

results = []

for sample in samples:

    sj = f"/mnt/faruq2/lab_users/shahrumi/RNA-Seq/SCA12_DIV30_neurons_rna_seq/RNA_Fastq_Dragen_Mapping/{sample}.SJ.out.tab"

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

out = "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/results/PPP2R2B_junctions_by_sample.tsv"

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
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/results/PPP2R2B_junctions_by_sample.tsv",
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
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/results/PPP2R2B_junction_summary.tsv",
    sep="\t", index=False
)

# Exploratory priority set, not a definitive validation filter.
candidates = summary[
    summary["novel_in_both"]
    & (summary["total_unique_reads"] >= 5)
    & (summary["max_overhang"] >= 10)
]

candidates.to_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/day30/results/PPP2R2B_prioritized_novel_junctions.tsv",
    sep="\t", index=False
)

print(candidates.to_string(index=False))
