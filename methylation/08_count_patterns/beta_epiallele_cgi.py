#!/usr/bin/env python3

import sys
import pandas as pd

infile, outfile = sys.argv[1], sys.argv[2]

df = pd.read_csv(infile, sep="\t", low_memory=False)

# CpG island identifier from island coordinates
df["cpg_island_id"] = (
    df["chrom"].astype(str) + ":"
    + df["start_region"].astype(str) + "-"
    + df["end_region"].astype(str)
)

# Force numeric columns and drop rows where conversion failed
df["nbase"] = pd.to_numeric(df["nbase"], errors="coerce")
df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
df = df.dropna(subset=["nbase", "beta"])

# Collapse from one row per read x CpG position to one row per read
reads = df[
    ["sample", "cpg_island_id", "pattern",
     "seqnames", "start", "end", "strand", "nbase", "beta"]
].drop_duplicates()

# Reads per sample x island x pattern
counts = (
    reads
    .groupby(["sample", "cpg_island_id", "pattern"])
    .agg(
        count=("pattern", "size"),
        mean_nbase=("nbase", "mean"),
        mean_beta=("beta", "mean"),
        sd_beta=("beta", "std"),
    )
    .reset_index()
)

# Total reads per sample x island
counts["total_reads"] = (
    counts.groupby(["sample", "cpg_island_id"])["count"].transform("sum")
)

counts.to_csv(outfile, sep="\t", index=False)
