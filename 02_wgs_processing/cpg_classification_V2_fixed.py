#!/usr/bin/env python3

import argparse
import pysam
from intervaltree import IntervalTree, Interval
from collections import defaultdict
from tqdm import tqdm

# ---------------------------
# Load CpG flanking sites
# ---------------------------
def load_flanking_sites(bed_path):
    trees = defaultdict(IntervalTree)
    total = 0
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            region, seq = line.split("\t")
            contig, coords = region.split(":")
            if contig.startswith("PUMA_"):
                contig = contig.replace("PUMA_", "", 1)
            start, end = map(int, coords.split("-"))
            trees[contig].add(Interval(start, end + 1, seq))
            total += 1
    print(f"Loaded {total} CpG flanking sites across {len(trees)} contigs", flush=True)
    return trees

# ---------------------------
# CpG classification
# ---------------------------
def classify_cpg(flank_seq, ref, alt):
    s1, s2, s3 = flank_seq[0], flank_seq[1], flank_seq[2]
    if s2 == "C" and s3 == "G":
        return "cpg1"
    elif s1 == "C" and s2 == "G":
        return "cpg2"
    elif alt == "C" and s3 == "G":
        return "poly-cpg3"
    elif alt == "G" and s1 == "C":
        return "poly-cpg4"
    else:
        return "variant"

# ---------------------------
# Process a single variant
# ---------------------------
def process_variant(rec, trees, samples, out_handle, cpg_summary):
    contig = rec.contig
    pos = rec.pos  # 1-based
    if contig not in trees:
        return
    overlaps = trees[contig].at(pos)
    if not overlaps:
        return

    alt_allele = rec.alts[0] if rec.alts else "."

    for interval in overlaps:
        status = classify_cpg(interval.data, rec.ref, alt_allele)

        row = {
            "contig": contig,
            "pos": pos,
            "ref": rec.ref,
            "alt": alt_allele,
            "cpg_seq": interval.data,
            "status": status,
        }

        for s in samples:
            gt = rec.samples[s]["GT"]
            if gt is None:
                # completely missing genotype
                row[s] = "NA"
                cpg_summary[(status, s)]["NA"] += 1
            else:
                # sum alleles only if they exist
                allele_sum = sum(a for a in gt[:2] if a is not None)
                allele_sum = min(allele_sum, 2)
                row[s] = allele_sum
                cpg_summary[(status, s)][allele_sum] += 1

        # write row immediately
        out_handle.write(
            "\t".join(str(row[h]) for h in ["contig","pos","ref","alt","cpg_seq","status"] + samples) + "\n"
        )

# ---------------------------
# Main
# ---------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--flanking", required=True)
    parser.add_argument("--out", required=True, help="Raw genotype TSV output")
    parser.add_argument("--summary", required=True, help="CpG allele count summary per individual")
    args = parser.parse_args()

    print("ARGS:", vars(args), flush=True)

    trees = load_flanking_sites(args.flanking)
    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    print(f"Processing VCF: {args.vcf} with {len(samples)} samples", flush=True)

    with open(args.out, "w") as out_handle:
        header = ["contig", "pos", "ref", "alt", "cpg_seq", "status"] + samples
        out_handle.write("\t".join(header) + "\n")

        cpg_summary = defaultdict(lambda: defaultdict(int))

        for rec in tqdm(vcf, desc="Reading VCF", unit="variants"):
            process_variant(rec, trees, samples, out_handle, cpg_summary)

    print("Raw genotype output complete.", flush=True)

    with open(args.summary, "w") as sum_out:
        sum_header = ["cpg_status", "sample", "allele", "count"]
        sum_out.write("\t".join(sum_header) + "\n")
        for (status, sample), allele_counts in cpg_summary.items():
            for allele, count in allele_counts.items():
                sum_out.write(f"{status}\t{sample}\t{allele}\t{count}\n")

    print("CpG summary output complete.", flush=True)


if __name__ == "__main__":
    main()

