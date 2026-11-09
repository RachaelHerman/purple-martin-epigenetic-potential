#!/usr/bin/env python3
import os
import sys
import multiprocessing as mp
from subprocess import Popen

SPLIT_FASTQ_SCRIPT = "split_fastq.py"

if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <file_list> <outdirectory> <nb_reads> <nbthreads>")
    sys.exit(1)

file_list = sys.argv[1]
outdirectory = sys.argv[2]
nb_reads = int(sys.argv[3])
nbthreads = int(sys.argv[4])

fileout_path = os.path.join(outdirectory, os.path.basename(file_list) + '.FINISHED')
fileout = open(fileout_path, 'w')

with open(file_list, 'r') as f:
    data = [line.strip() for line in f if line.strip()]

runs = []
samp_dirs = {}

for filepath in data:
    if filepath.endswith('.R1.fastq.gz'):
        parts = filepath.split('/')
        directory = '/'.join(parts[:-1])
        samp_dir = parts[-2]
        run_id = parts[-1].replace('.R1.fastq.gz', '')
        if samp_dir not in samp_dirs:
            samp_dirs[samp_dir] = []
        samp_dirs[samp_dir].append(run_id)
        batch_name = parts[-2].split('_')[1] + '_' + str(samp_dirs[samp_dir].index(run_id) + 1)
        runs.append([directory, samp_dir, run_id, batch_name])

    elif filepath.endswith('_1.fq.gz'):
        parts = filepath.split('/')
        directory = '/'.join(parts[:-1])
        samp_dir = parts[-2]
        run_id = parts[-1].replace('_1.fq.gz', '')
        if samp_dir not in samp_dirs:
            samp_dirs[samp_dir] = []
        samp_dirs[samp_dir].append(run_id)
        batch_name = samp_dir + '_' + str(samp_dirs[samp_dir].index(run_id) + 1)
        runs.append([directory, samp_dir, run_id, batch_name])

for samp_dir in samp_dirs.keys():
    os.makedirs(os.path.join(outdirectory, samp_dir), exist_ok=True)

batches = [runs[i:i + nbthreads] for i in range(0, len(runs), nbthreads)]
cwd = os.getcwd()

def split_fastq(idx, splits, outdirectory, nb_reads, output_queue):
    os.chdir(cwd)
    directory, samp_dir, run_id, batch_name = splits[idx]
    output_path = os.path.join(outdirectory, samp_dir, batch_name)
    cmd = f"{SPLIT_FASTQ_SCRIPT} {run_id} {directory} {output_path} {nb_reads}"
    proc = Popen(cmd, shell=True)
    proc.wait()
    output_queue.put(f"finished {batch_name}")

for batch in batches:
    output_queue = mp.Queue()
    processes = [mp.Process(target=split_fastq, args=(i, batch, outdirectory, nb_reads, output_queue))
                 for i in range(len(batch))]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    results = [output_queue.get() for _ in processes]
    for r in results:
        fileout.write(r + '\n')

fileout.close()
