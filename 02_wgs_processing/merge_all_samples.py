#!/usr/bin/env python3
import os
import re
import argparse
import subprocess
from glob import glob

parser = argparse.ArgumentParser(description="Merge per-contig GVCFs for multiple samples")
parser.add_argument("--indir", required=True, help="Directory containing intermediate gVCFs")
parser.add_argument("--outdir", required=True, help="Output directory for final merged gVCFs")
parser.add_argument("--gatk", required=True, help="Path to GATK jar")
parser.add_argument("--java", required=True, help="Path to Java binary")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# Regex for filenames like:
# PM223_PE_ME_merged_JAKJHU010002895.1_JAKJHU010000029.1.g.vcf.gz
pattern = re.compile(r"^(?P<sample>.+?)_JAKJHU.*\.g\.vcf\.gz$")

files = glob(os.path.join(args.indir, "*.g.vcf.gz"))

# Group files by sample
samples = {}
for f in files:
    base = os.path.basename(f)
    m = pattern.match(base)
    if not m:
        continue
    s = m.group("sample")
    samples.setdefault(s, []).append(f)

if not samples:
    raise SystemExit("ERROR: No matching gVCF files found.")

for sample, vcfs in samples.items():
    print(f"\n=== Processing sample: {sample} ===")
    log_path = os.path.join(args.outdir, f"{sample}_merge_log.txt")

    # Sort by contig numerically where possible
    def contig_sort_key(path):
        m = re.search(r"_(JAKJHU[\w\.]+)_(JAKJHU[\w\.]+)\.g\.vcf\.gz$", path)
        return (m.group(1), m.group(2)) if m else path

    vcfs_sorted = sorted(vcfs, key=contig_sort_key)

    # Write merge log
    with open(log_path, "w") as log:
        log.write(f"Sample: {sample}\n")
        log.write(f"Total contig batches: {len(vcfs_sorted)}\n\n")
        log.write("Included files:\n")
        for v in vcfs_sorted:
            log.write(os.path.basename(v) + "\n")

    final_vcf = os.path.join(args.outdir, f"{sample}.g.vcf.gz")
    print(f"Merging into: {final_vcf}")

    merge_cmd = [
        args.java, "-Xmx60G", "-jar", args.gatk,
        "GatherVcfs"
    ]
    for v in vcfs_sorted:
        merge_cmd.extend(["-I", v])
    merge_cmd.extend(["-O", final_vcf])

    print("Running merge command...")
    subprocess.run(merge_cmd, check=True)

    # Index final VCF
    print("Indexing final VCF...")
    subprocess.run(["gatk", "IndexFeatureFile", "-I", final_vcf], check=True)

    # Delete intermediate files
    print("Cleaning up intermediate files...")
    for v in vcfs_sorted:
        try:
            os.remove(v)
        except:
            pass
        if os.path.exists(v + ".tbi"):
            os.remove(v + ".tbi")

    print(f"✔ Done: {sample}")
    print(f"✔ Log written to: {log_path}\n")

