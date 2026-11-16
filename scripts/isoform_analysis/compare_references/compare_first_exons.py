python - <<'PY'
import pandas as pd

files = [
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.gencode.catalog.tsv",
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.refseq.catalog.tsv"
]

df = pd.concat(
    [pd.read_csv(f, sep="\t") for f in files],
    ignore_index=True
)

summary = (
    df.groupby(["chromosome", "strand", "first_exon"])
    .agg(
        transcript_count=("transcript_id", "nunique"),
        transcript_ids=("transcript_id",
                        lambda x: ",".join(map(str, x))),
        sources=("source",
                 lambda x: ",".join(sorted(set(x))))
    )
    .reset_index()
)

summary.to_csv(
    "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B_first_exons.tsv",
    sep="\t",
    index=False
)

print(summary.to_string(index=False))
PY
