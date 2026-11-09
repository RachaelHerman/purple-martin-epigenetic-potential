#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
from multiprocessing import Pool
import argparse
import tempfile

parser = argparse.ArgumentParser(description="Parallelized GATK HaplotypeCaller on a single BAM")
parser.add_argument("--bam", required=True)
parser.add_argument("--reference", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--threads", type=int, default=160)
parser.add_argument("--batch_size", type=int, default=20)
parser.add_argument("--threads_per_batch", type=int, default=4)
parser.add_argument("--gatk", required=True)
parser.add_argument("--java", required=True)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
sample = os.path.basename(args.bam).split('.')[0]

contigs = [r.id for r in SeqIO.parse(args.reference, "fasta")]
batches = [contigs[i:i+args.batch_size] for i in range(0, len(contigs), args.batch_size)]

def run_batch(batch):
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".list", dir=args.outdir) as f:
        for contig in batch:
            f.write(contig + "\n")
        interval_file = f.name

    out_vcf = os.path.join(args.outdir, f"{sample}_{batch[0]}_{batch[-1]}.g.vcf.gz")
    cmd = [
        args.java, "-Xmx60G", "-jar", args.gatk,
        "HaplotypeCaller",
        "-R", args.reference,
        "-I", args.bam,
        "--emit-ref-confidence", "GVCF",
        "-L", interval_file,
        "-O", out_vcf,
        "--native-pair-hmm-threads", str(args.threads_per_batch)
    ]
    subprocess.run(cmd, check=True)
    os.remove(interval_file)
    return out_vcf

max_parallel_batches = args.threads // args.threads_per_batch
with Pool(processes=max_parallel_batches) as pool:
    vcf_files = pool.map(run_batch, batches)

merged_vcf = os.path.join(args.outdir, f"{sample}.g.vcf.gz")
merge_cmd = [
    args.java, "-Xmx60G", "-jar", args.gatk,
    "GatherVcfs"
] + [f"-I {vcf}" for vcf in vcf_files] + ["-O", merged_vcf]
subprocess.run(" ".join(merge_cmd), shell=True, check=True)

for vcf in vcf_files:
    os.remove(vcf)
    if os.path.exists(vcf + ".tbi"):
        os.remove(vcf + ".tbi")
