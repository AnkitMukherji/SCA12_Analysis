#!/usr/bin/env python3

# EXACT_ANNOTATED: matches a known annotated exon chain.
# KNOWN_INTRON_CHAIN_DIFFERENT_TERMINI: same splicing pattern, but differing transcript boundaries.
# NOVEL_SPLICE_JUNCTION: at least one splice junction absent from both reference catalogues.
# NEW_COMBINATION_OF_KNOWN_JUNCTIONS: a new intron-chain combination constructed from already annotated junctions.
# SINGLE_EXON_UNRESOLVED: cannot establish a multi-exon transcript structure.

import re
import csv
import os
from collections import defaultdict

# Reference Catalogs (TSVs)
GENCODE_CATALOG = "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.gencode.catalog.tsv"
REFSEQ_CATALOG = "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/compare_ncbi_gencode_transcripts/PPP2R2B.refseq.catalog.tsv"

# The two merged assemblies generated from StringTie
ASSEMBLIES = {
    "gencode_guided": "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/annotation/PPP2R2B_assembled_gencode.gtf",
    "refseq_guided": "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/annotation/PPP2R2B_assembled_refseq.gtf"
}

def read_catalog(filename):
    if not os.path.exists(filename):
        return []
    with open(filename) as f:
        return list(csv.DictReader(f, delimiter="\t"))

# Merge catalogs for cross-reference matching
refs = read_catalog(GENCODE_CATALOG) + read_catalog(REFSEQ_CATALOG)

def attrs(s):
    return dict(re.findall(r'([\w.:-]+)\s+"([^"]*)"', s))

# Load reference PPP2R2B exons to avoid assigning all
# transcripts overlapping this large genomic locus to the gene.
ref_exons = []
for filename in ["/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/annotation/PPP2R2B_assembled_gencode.gtf", "/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/annotation/PPP2R2B_assembled_refseq.gtf"]:
    if os.path.exists(filename):
        with open(filename) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                v = line.rstrip("\n").split("\t")
                if len(v) == 9 and v[2] == "exon":
                    ref_exons.append((v[0], v[6], int(v[3]), int(v[4])))

def classify_assembly(assembly_path, output_path):
    # Collect StringTie exon structures.
    models = defaultdict(list)
    model_info = {}

    if not os.path.exists(assembly_path):
        print(f"Warning: Assembly file not found: {assembly_path}")
        return

    with open(assembly_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            v = line.rstrip("\n").split("\t")
            if len(v) != 9:
                continue

            a = attrs(v[8])
            tid = a.get("transcript_id")
            if not tid:
                continue

            model_info[tid] = (v[0], v[6], a.get("gene_id", ""))

            if v[2] == "exon":
                models[tid].append((int(v[3]), int(v[4])))

    rows = []

    for tid, raw_exons in sorted(models.items()):
        chrom, strand, gene_id = model_info[tid]
        exons = sorted(set(raw_exons))

        # This primary analysis concerns the PPP2R2B sense strand (-).
        if strand != "-":
            continue

        # Require actual exon overlap with annotated PPP2R2B exons.
        shared = any(
            chrom == rc
            and strand == rs
            and s1 <= e2
            and s2 <= e1
            for s1, e1 in exons
            for rc, rs, s2, e2 in ref_exons
        )

        if not shared:
            continue

        introns = [
            (exons[i][1] + 1, exons[i+1][0] - 1)
            for i in range(len(exons)-1)
        ]

        fmt = lambda values: ";".join(f"{s}-{e}" for s, e in values)

        exon_chain = fmt(exons)
        intron_chain = fmt(introns)
        first_exon = exons[-1]   # Negative strand (PPP2R2B)

        compatible_refs = [
            r for r in refs
            if r["chromosome"] == chrom
            and r["strand"] == strand
        ]

        exact = [r for r in compatible_refs if r["exons"] == exon_chain]
        same_introns = [r for r in compatible_refs if introns and r["introns"] == intron_chain]

        known_junctions = set()
        for r in compatible_refs:
            if r["introns"]:
                known_junctions.update(r["introns"].split(";"))

        novel_junctions = [
            j for j in introns
            if fmt([j]) not in known_junctions
        ]

        if exact:
            category = "EXACT_ANNOTATED"
        elif same_introns:
            category = "KNOWN_INTRON_CHAIN_DIFFERENT_TERMINI"
        elif novel_junctions:
            category = "NOVEL_SPLICE_JUNCTION"
        elif introns:
            category = "NEW_COMBINATION_OF_KNOWN_JUNCTIONS"
        else:
            category = "SINGLE_EXON_UNRESOLVED"

        matches = exact if exact else same_introns

        rows.append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "category": category,
            "reference_matches": ",".join(
                f'{r["source"]}:{r["transcript_id"]}'
                for r in matches
            ),
            "first_exon": fmt([first_exon]),
            "exon_count": len(exons),
            "exons": exon_chain,
            "introns": intron_chain,
            "novel_junctions": fmt(novel_junctions)
        })

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="") as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    print(f"Classified {len(rows)} sense-strand models from {os.path.basename(assembly_path)}")
    print(f"Saved to: {output_path}\n")


# Run classification for both assemblies
for name, path in ASSEMBLIES.items():
    out_tsv = f"/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_short_read/results/PPP2R2B_{name}_classification.tsv"
    classify_assembly(path, out_tsv)
