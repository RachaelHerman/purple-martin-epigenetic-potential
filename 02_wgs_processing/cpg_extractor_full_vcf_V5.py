#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF
import pandas as pd
from intervaltree import IntervalTree
from collections import defaultdict
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Extract CpG genotypes from full VCF")
    parser.add_argument("--vcf", required=True, help="Full multi-contig VCF file")
    parser.add_argument("--flanking", required=True, help="CpG flanking sites BED-like file")
    parser.add_argument("--chunk_size", type=int, default=10_000_000, help="Number of variants per chunk")
    parser.add_argument("--out", required=True, help="Output file path")
    return parser.parse_args()

def load_flanking_sites(flanking_path):
    print(f"Loading CpG flanking sites from {flanking_path}...")
    df = pd.read_csv(flanking_path, sep="\t", header=None, names=["region", "sequence"], dtype=str)
    contig_coords = df["region"].str.split("[:-]", expand=True)
    df["contig"] = contig_coords[0]
    df["start"] = contig_coords[1].astype(int)
    df["end"] = contig_coords[2].astype(int)

    trees = defaultdict(IntervalTree)
    for contig, group in df.groupby("contig"):
        intervals = [(row.start, row.end, row.sequence) for row in group.itertuples(index=False)]
        trees[contig].update(intervals)

    print(f"Built interval trees for {len(trees)} contigs with {len(df)} CpG sites")
    return trees

def process_vcf(vcf_path, trees, chunk_size, out_path):
    vcf = VCF(vcf_path)
    samples = vcf.samples
    print(f"Processing VCF: {vcf_path} in chunks of {chunk_size} variants")
    print(f"Found {len(samples)} samples")

    with open(out_path, "w") as fout:
        fout.write("contig\tpos\tref\talt\t" + "\t".join(samples) + "\n")
        chunk = []
        for i, var in enumerate(vcf):
            chunk.append(var)
            if len(chunk) >= chunk_size:
                write_chunk(chunk, trees, fout)
                chunk = []

        if chunk:
            write_chunk(chunk, trees, fout)

def write_chunk(chunk, trees, fout):
    for var in chunk:
        contig = var.CHROM
        pos = var.POS
        ref = var.REF
        alts = var.ALT
        # Only process if contig has CpG sites
        if contig not in trees:
            continue
        overlaps = trees[contig].overlap(pos, pos+1)
        if not overlaps:
            continue
        # Write a row per variant with genotypes for all samples
        gts = [var.genotypes[j][:2] for j in range(len(var.genotypes))]  # [GT1, GT2]
        gt_strs = ["/".join(map(str, gt)) for gt in gts]
        fout.write(f"{contig}\t{pos}\t{ref}\t{','.join(alts)}\t" + "\t".join(gt_strs) + "\n")

if __name__ == "__main__":
    args = parse_args()
    print(f"ARGS: {vars(args)}")

    flanking_trees = load_flanking_sites(args.flanking)
    process_vcf(args.vcf, flanking_trees, args.chunk_size, args.out)
    print(f"Finished writing output to {args.out}")

