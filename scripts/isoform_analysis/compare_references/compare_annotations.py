#!/usr/bin/env python3

import pandas as pd

g = pd.read_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.gencode.catalog.tsv",
    sep="\t", dtype=str
).fillna("")

r = pd.read_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.refseq.catalog.tsv",
    sep="\t", dtype=str
).fillna("")

rows = []

for _, a in g.iterrows(): # Skipping index using _
    for _, b in r.iterrows(): # Skipping index using _

        if (
            a["chromosome"] != b["chromosome"] or a["strand"] != b["strand"]
        ):
            continue

        if a["exons"] == b["exons"]:
            relationship = "EXACT_EXON_CHAIN"

        elif (
            a["introns"] and a["introns"] == b["introns"]
        ):
            relationship = "SAME_INTRON_CHAIN_DIFFERENT_ENDS"

        elif (
            a["introns"] and b["introns"] and set(a["introns"].split(";")) & set(b["introns"].split(";"))
        ):
            relationship = "PARTIAL_INTRON_OVERLAP"

        else:
            continue

        rows.append({
            "GENCODE": a["transcript_id"],
            "REFSEQ": b["transcript_id"],
            "relationship": relationship,
            "GENCODE_first_exon": a["first_exon"],
            "REFSEQ_first_exon": b["first_exon"]
        })

pd.DataFrame(
    rows,
    columns=[
        "GENCODE", "REFSEQ", "relationship",
        "GENCODE_first_exon", "REFSEQ_first_exon"
    ]
).to_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/GENCODE_REFSEQ_matches.tsv",
    sep="\t", index=False
)

print("Compared GENCODE and RefSeq transcript structures.")
