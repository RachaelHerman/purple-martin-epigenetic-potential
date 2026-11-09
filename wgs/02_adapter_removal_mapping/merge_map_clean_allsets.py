#!/usr/bin/env python3
import os
import sys
from subprocess import Popen
import multiprocessing as mp
import re

def find_sample_sets(sample_name, sample_dir):
    files = os.listdir(sample_dir)
    pattern = re.compile(rf"^{sample_name}_(\d+)_split_\d+_R1\.fq\.gz$")
    sets = set()
    for f in files:
        match = pattern.match(f)
        if match:
            sets.add(match.group(1))
    return sorted(list(sets), key=int)

def find_split_files(sample_name, read_set, sample_dir):
    files = os.listdir(sample_dir)
    split_files = []
    for f in files:
        if f.startswith(f"{sample_name}_{read_set}_split") and f.endswith('R1.fq.gz'):
            parts = f.split('_')
            try:
                split_files.append((int(parts[3]), f))
            except ValueError:
                continue
    split_files.sort()
    return split_files

def process_split(x, batch, sample_name, read_set, sample_dir, bwa_index_prefix, nbthreads, bams_dir, output_queue):
    split_nb, _ = batch[x]
    targ_stem = f"{sample_name}_{read_set}_split_{split_nb}"

    adaptrem_dir = os.path.join(sample_dir, 'adaptrem')
    os.makedirs(adaptrem_dir, exist_ok=True)
    os.makedirs(bams_dir, exist_ok=True)

    r1 = os.path.join(sample_dir, f"{targ_stem}_R1.fq.gz")
    r2 = os.path.join(sample_dir, f"{targ_stem}_R2.fq.gz")
    basename = os.path.join(adaptrem_dir, targ_stem)

    cmd_adapter = (
        f"AdapterRemoval --file1 {r1} --file2 {r2} "
        f"--basename {basename} --trimns --trimqualities --collapse --gzip"
    )
    Popen(cmd_adapter, shell=True).wait()

    cmds_bwa = [
        f"bwa mem -M -t {nbthreads} {bwa_index_prefix} {basename}.pair1.truncated.gz {basename}.pair2.truncated.gz | samtools view -Sb - > {bams_dir}/{targ_stem}_PE.bam",
        f"bwa mem -M -t {nbthreads} {bwa_index_prefix} {basename}.collapsed.gz | samtools view -Sb - > {bams_dir}/{targ_stem}_ME.bam",
        f"bwa mem -M -t {nbthreads} {bwa_index_prefix} {basename}.collapsed.truncated.gz | samtools view -Sb - > {bams_dir}/{targ_stem}_MEt.bam"
    ]
    for cmd in cmds_bwa:
        Popen(cmd, shell=True).wait()

    output_queue.put(f"Finished {targ_stem}")

def main():
    sample_name = sys.argv[1]
    sample_dir = sys.argv[2]
    nbthreads = max(1, int(sys.argv[3]) // 2)
    bwa_index_prefix = sys.argv[4]

    bams_dir = os.path.join(sample_dir, "adaptrem", "bams")

    for read_set in find_sample_sets(sample_name, sample_dir):
        split_files = find_split_files(sample_name, read_set, sample_dir)
        batches = [split_files[i:i + nbthreads] for i in range(0, len(split_files), nbthreads)]
        for batch in batches:
            output_queue = mp.Queue()
            processes = [
                mp.Process(
                    target=process_split,
                    args=(i, batch, sample_name, read_set, sample_dir, bwa_index_prefix, nbthreads, bams_dir, output_queue)
                )
                for i in range(len(batch))
            ]
            for p in processes: p.start()
            for p in processes: p.join()
            [output_queue.get() for _ in processes]

if __name__ == "__main__":
    main()
