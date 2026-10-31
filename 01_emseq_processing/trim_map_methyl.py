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
nbthreads = nbthreads_init // 2
ref = argv[4]            # bismark reference directory

cwd = os.getcwd()
trim_env = "trim_env"
bismark_env = "bismark_env"

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
batches = [use_files[i:i + nbthreads] for i in range(0, nb_process, nbthreads)]

# -----------------------------
# Worker functions
# -----------------------------
def trim_worker(split):
    split_num, filename = split
    prefix = filename.split("_split")[0]
    targ_stem = f"{prefix}_split_{split_num}"
    trim_dir = os.path.join(cwd, samp_dir, "trim_galore")
    os.makedirs(trim_dir, exist_ok=True)

    trim_cmd = (
        f"trim_galore "
        f"--output_dir {trim_dir} "
        f"--paired --cores 4 --quality 20 --fastqc "
        f"--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        f"--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 "
        f"{targ_stem}_R1.fq.gz {targ_stem}_R2.fq.gz"
    )
    subprocess.run(trim_cmd, shell=True, check=True)
    return targ_stem

def bismark_worker(targ_stem):
    trim_dir = os.path.join(cwd, samp_dir, "trim_galore")
    bam_dir = os.path.join(trim_dir, "bam", "bams")
    os.makedirs(bam_dir, exist_ok=True)

    tr1 = os.path.join(trim_dir, f"{targ_stem}_R1_val_1.fq.gz")
    tr2 = os.path.join(trim_dir, f"{targ_stem}_R2_val_2.fq.gz")

    bismark_cmd = [
        "bismark",
        "-o", bam_dir,
        ref,
        "-1", tr1,
        "-2", tr2
    ]
    subprocess.run(bismark_cmd, check=True)
    return f"finished {targ_stem}"

# -----------------------------
# Run batches
# -----------------------------
for batch in batches:
    # Activate trim environment once for the batch
    trim_cmd_prefix = f"source ~/.bashrc && conda activate {trim_env} && "
    pool = mp.Pool(processes=len(batch))
    trim_results = pool.map(lambda s: subprocess.run(trim_cmd_prefix + f"trim_galore ...", shell=True), batch)
    pool.close()
    pool.join()

    # Activate bismark environment once for the batch
    bismark_cmd_prefix = f"source ~/.bashrc && conda activate {bismark_env} && "
    pool = mp.Pool(processes=len(batch))
    bismark_results = pool.map(lambda t: subprocess.run(bismark_cmd_prefix + f"bismark ...", shell=True), trim_results)
    pool.close()
    pool.join()

    print(bismark_results)

