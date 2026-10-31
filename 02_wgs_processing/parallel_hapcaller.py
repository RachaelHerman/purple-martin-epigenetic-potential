#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
from multiprocessing import Pool
import argparse
import tempfile

parser = argparse.ArgumentParser(description="Parallelized GATK HaplotypeCaller on a single BAM")
parser.add_argument("--bam", required=True, help="Input BAM file")
parser.add_argument("--reference", required=True, help="Reference genome FASTA")
parser.add_argument("--outdir", required=True, help="Output directory for GVCFs")
parser.add_argument("--threads", type=int, default=160, help="Total threads available")
parser.add_argument("--batch_size", type=int, default=20, help="Number of contigs per batch")
parser.add_argument("--threads_per_batch", type=int, default=4, help="Threads to use per HaplotypeCaller batch")
parser.add_argument("--gatk", required=True, help="Path to GATK jar")
parser.add_argument("--java", required=True, help="Path to Java binary")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
sample = os.path.basename(args.bam).split('.')[0]

# Get contigs
contigs = [r.id for r in SeqIO.parse(args.reference, "fasta")]

# Split contigs into batches
batches = [contigs[i:i+args.batch_size] for i in range(0, len(contigs), args.batch_size)]

def run_batch(batch):
    # Write temporary interval file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".list", dir=args.outdir) as f:
        for contig in batch:
            f.write(contig + "\n")
        interval_file = f.name

    out_vcf = os.path.join(args.outdir, f"{sample}_{batch[0]}_{batch[-1]}.g.vcf.gz")
    cmd = [
        args.java, f"-Xmx60G", "-jar", args.gatk,
        "HaplotypeCaller",
        "-R", args.reference,
        "-I", args.bam,
        "--emit-ref-confidence", "GVCF",
        "-L", interval_file,
        "-O", out_vcf,
        "--native-pair-hmm-threads", str(args.threads_per_batch)
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Remove temporary interval file
    os.remove(interval_file)

    return out_vcf

max_parallel_batches = args.threads // args.threads_per_batch
with Pool(processes=max_parallel_batches) as pool:
    vcf_files = pool.map(run_batch, batches)

# Merge all batch GVCFs
merged_vcf = os.path.join(args.outdir, f"{sample}.g.vcf.gz")
merge_cmd = [
    args.java, "-Xmx60G", "-jar", args.gatk,
    "GatherVcfs"
] + [f"-I {vcf}" for vcf in vcf_files] + ["-O", merged_vcf]
print("Merging GVCFs into:", merged_vcf)
subprocess.run(" ".join(merge_cmd), shell=True, check=True)

# Optional: cleanup intermediate batch GVCFs
for vcf in vcf_files:
    os.remove(vcf)
    if os.path.exists(vcf + ".tbi"):
        os.remove(vcf + ".tbi")

print(f"Final merged GVCF: {merged_vcf}")

