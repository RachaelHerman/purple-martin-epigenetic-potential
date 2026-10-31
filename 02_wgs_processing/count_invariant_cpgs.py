#!/usr/bin/env python3

import pandas as pd

# Load your output file
df = pd.read_csv("PUMA_allsites_cpg_class_genotypes.tsv", sep="\t")

# Consider all CpG types
cpg_types = ['cpg1', 'cpg2', 'polycpg3', 'polycpg4']

cpg_df = df[df['status'].isin(cpg_types)]

# Identify sample columns (adjust if needed)
sample_cols = df.columns[6:]  # assuming first 6 cols are contig,pos,ref,alt,cpg_seq,status

# Filter invariant CpGs (all individuals homozygous ref or 0)
# This works even if polycpgs are present
invariant_cpgs = cpg_df[(cpg_df[sample_cols] == 0).all(axis=1)]

print(f"Found {len(invariant_cpgs)} invariant CpGs (including polyCpGs)")
print(invariant_cpgs[['contig','pos','ref','alt','cpg_seq','status']].head())

# Optionally, write to a new file
invariant_cpgs.to_csv("invariant_cpgs.tsv", sep="\t", index=False)

