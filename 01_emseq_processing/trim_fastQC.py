#!/usr/bin/env python3
import os
import subprocess
import multiprocessing as mp
from sys import argv

# -----------------------------
# Arguments
# -----------------------------
samp_name = argv[1]      # sample prefix, e.g. AZH055_E
samp_dir = argv[2]       # directory containing split reads
nbthreads_init = int(argv[3])
nbthreads = nbthreads_init // 2  # threads per batch

cwd = os.getcwd()
trim_env = "trim_env"    # conda environment with trim_galore

# -----------------------------
# Gather split files
# -----------------------------
files = os.listdir(os.path.join(cwd, samp_dir))
use_files = []

for f in files:
    if f.startswith(samp_name) and f.endswith("_R1.fq.gz") and "_split" in f:
        try:
            split_num = int(f.split("_split")[1].split("_")[0])
            use_files.append([split_num, f])
        except ValueError:
            pass

use_files.sort()
nb_process = len(use_files)

# Split into batches for multiprocessing
batches = [use_files[i:i + nbthreads] for i in range(0, nb_process, nbthreads)]

# -----------------------------
# Worker function
# -----------------------------
def trim_worker(split):
    split_num, filename = split
    prefix = filename.split("_split")[0]
    targ_stem = f"{prefix}_split_{split_num}"

    trim_dir = os.path.join(cwd, samp_dir, "trim_galore")
    os.makedirs(trim_dir, exist_ok=True)

    trim_cmd = (
        f"source ~/.bashrc && conda activate {trim_env} && "
        f"trim_galore "
        f"--output_dir {trim_dir} "
        f"--paired --cores 4 --quality 20 --fastqc "
        f"--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        f"--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 "
        f"{targ_stem}_R1.fq.gz {targ_stem}_R2.fq.gz"
    )
    print(f"Running: {trim_cmd}")
    subprocess.run(trim_cmd, shell=True, check=True)
    return f"finished {targ_stem}"

# -----------------------------
# Run batches
# -----------------------------
for batch in batches:
    pool = mp.Pool(processes=len(batch))
    results = pool.map(trim_worker, batch)
    pool.close()
    pool.join()
    print(results)

