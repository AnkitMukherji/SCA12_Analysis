
#!/usr/bin/env python3

import re
import csv
import sys
from collections import defaultdict

# Run with:
# python reference_catalog.py input.gtf source output_prefix

input_gtf, source, prefix = sys.argv[1:4]

TARGET = "PPP2R2B"
GENE_ID = {"5521", "GeneID:5521"}

def attributes(text):
    """Parse standard GTF key-value attributes."""
    return dict(re.findall(r'([\w.:-]+)\s+"([^"]*)"', text))

def is_target(a):
    return (
        a.get("gene_name") == TARGET
        or a.get("gene") == TARGET
        or a.get("gene_id") in GENE_ID
        or "GeneID:5521" in a.get("db_xref", "").split(",")
    )

records = []
target_ids = set()

# Extract transcript IDs associated with PPP2R2B.
with open(input_gtf) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        a = attributes(fields[8])
        tid = a.get("transcript_id")
        records.append((line, fields, a))
        if tid and is_target(a):
            target_ids.add(tid)

if not target_ids:
    sys.exit(f"No PPP2R2B transcripts found in {input_gtf}")

# Retain all features associated with identified transcripts.
exons = defaultdict(list)
cds = defaultdict(list)
strands = {}
chromosomes = {}

with open(prefix + ".gtf", "w") as out:
    for line, fields, a in records:
        tid = a.get("transcript_id")
        if tid not in target_ids:
            continue

        out.write(line)
        chrom, _, feature, start, end, _, strand, _, _ = fields

        chromosomes[tid] = chrom
        strands[tid] = strand

        if feature == "exon":
            exons[tid].append((int(start), int(end)))
        elif feature == "CDS":
            cds[tid].append((int(start), int(end)))

columns = [
    "source", "transcript_id", "chromosome", "strand",
    "start", "end", "exon_count", "first_exon",
    "first_intron", "exons", "introns", "cds"
]

with open(prefix + ".catalog.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(columns)

    for tid in sorted(target_ids):
        ee = sorted(set(exons[tid]))

        if not ee:
            continue

        # GTF coordinates are 1-based, inclusive.
        introns = [
            (ee[i][1] + 1, ee[i + 1][0] - 1)
            for i in range(len(ee) - 1)
        ]

        # For a minus-strand gene, the first exon is the
        # exon with the highest genomic coordinate.
        minus = strands[tid] == "-"
        first_exon = ee[-1] if minus else ee[0]
        first_intron = (
            introns[-1] if minus else introns[0]
        ) if introns else None

        fmt = lambda items: ";".join(
            f"{s}-{e}" for s, e in items
        )

        writer.writerow([
            source, tid, chromosomes[tid], strands[tid],
            ee[0][0], ee[-1][1],
            len(ee),
            fmt([first_exon]),
            fmt([first_intron]) if first_intron else "",
            fmt(ee),
            fmt(introns),
            fmt(sorted(set(cds[tid])))
        ])

print(f"{source}: {len(exons)} transcripts with exon records")
print(f"Created {prefix}.gtf")
print(f"Created {prefix}.catalog.tsv")
