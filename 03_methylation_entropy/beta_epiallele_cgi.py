#!/usr/bin/env python3

import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

print(f"Processing {infile}")

df = pd.read_csv(
    infile,
    sep="\t",
    low_memory=False
)

# make CpG island identifier
df["promoter_id"] = (
    df["chrom"].astype(str)
    + ":"
    + df["start_region"].astype(str)
    + "-"
    + df["end_region"].astype(str)
)

# force numeric columns
df["nbase"] = pd.to_numeric(df["nbase"], errors="coerce")
df["beta"] = pd.to_numeric(df["beta"], errors="coerce")

# remove rows where numeric values failed
df = df.dropna(subset=["nbase","beta"])

reads = df[
    [
        "sample",
        "promoter_id",
        "pattern",
        "seqnames",
        "start",
        "end",
        "strand",
        "nbase",
        "beta"
    ]
].drop_duplicates()

print(f"Unique reads retained: {len(reads)}")

counts = (
    reads
    .groupby(
        [
            "sample",
            "promoter_id",
            "pattern"
        ]
    )
    .agg(
        count=("pattern","size"),
        mean_nbase=("nbase","mean"),
        mean_beta=("beta","mean"),
        sd_beta=("beta","std")
    )
    .reset_index()
)

counts["total_reads"] = (
    counts
    .groupby(["sample","promoter_id"])
    ["count"]
    .transform("sum")
)

counts.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(f"Done: {outfile}")
