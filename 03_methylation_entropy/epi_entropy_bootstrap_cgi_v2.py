#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np


# Input Files

files = sys.argv[1].split()
depth_file = sys.argv[2]
out_file = sys.argv[3]
n_boot = int(sys.argv[4])

print(f"[INFO] files: {len(files)}", flush=True)


# Population assignment

north = {
    "AZH055", "AZH091", "FG05K08", "FG05K10", "FG07K03",
    "FG07K23", "FG08K10", "MB424", "MBA267", "ONB072"
}

south = {
    "FE13K09", "FE13K10", "FE13K36", "FE14K44", "FE14K57",
    "FE14K78", "FE14K96", "FE15K28", "FE15K33", "FE15K50"
}

ALL_SAMPLES = north | south


def extract_sample(filepath):

    basename = os.path.basename(filepath)
    for s in ALL_SAMPLES:
        if basename.startswith(s):
            return s
    return None


def get_pop(sample):
    if sample in north:
        return "north"
    if sample in south:
        return "south"
    return None


# Load CGI Depth Table

depth = pd.read_csv(
    depth_file,
    sep="\t"
)

required = [
    "cpg_island_id",
    "D_CpG_island"
]

for c in required:

    if c not in depth.columns:

        raise RuntimeError(
            f"Missing {c}. Found: {depth.columns.tolist()}"
        )


depth["D_CpG_island"] = pd.to_numeric(
    depth["D_CpG_island"],
    errors="coerce"
)

depth = depth.dropna(
    subset=["D_CpG_island"]
)

depth = depth[
    depth["D_CpG_island"] >= 10
]

depth_map = dict(
    zip(
        depth["cpg_island_id"],
        depth["D_CpG_island"].astype(int)
    )
)

print(
    f"[INFO] usable CpG islands: {len(depth_map)}",
    flush=True
)


# Define entropy function

def entropy(patterns):

    _, counts = np.unique(
        patterns,
        return_counts=True
    )

    p = counts / counts.sum()

    return -np.sum(
        p * np.log2(p)
    )


# Load and compress data

tables = []


for f in files:

    sample = extract_sample(f)

    if sample is None:
        print(f"[WARNING] could not match {os.path.basename(f)} to a known sample -- skipping", flush=True)
        continue

    pop = get_pop(sample)

    if pop is None:
        continue

    print(
        f"[LOAD] {sample}",
        flush=True
    )

    df = pd.read_csv(
        f,
        sep="\t",
        usecols=[
            "promoter_id",
            "pattern"
        ],
        low_memory=False
    )

    df = df.rename(columns={"promoter_id": "cpg_island_id"})

    df = df[
        df["cpg_island_id"].isin(depth_map)
    ]

    if df.empty:
        continue

    counts = (
        df.groupby(
            [
                "cpg_island_id",
                "pattern"
            ]
        )
        .size()
        .reset_index(
            name="count"
        )
    )

    counts["sample"] = sample
    counts["population"] = pop

    tables.append(counts)

    del df
    del counts


read_table = pd.concat(
    tables,
    ignore_index=True
)

print(
    f"[INFO] grouped rows: {len(read_table)}",
    flush=True
)

print(
    f"[INFO] samples loaded: {sorted(read_table['sample'].unique())}",
    flush=True
)

del tables



# Bootstrapping regime (parallelized across islands)

import os as _os
from multiprocessing import Pool

workers = int(_os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"[INFO] parallel workers: {workers}", flush=True)

## ---- Pre-split read_table by island ONCE, shared to workers via initializer ----
ISLAND_TABLES = {pid: sub for pid, sub in read_table.groupby("cpg_island_id")}
DEPTH_MAP_G = None
NBOOT_G = None

def init_worker(island_tables, depth_map_arg, nboot_arg):
    global ISLAND_TABLES, DEPTH_MAP_G, NBOOT_G
    ISLAND_TABLES = island_tables
    DEPTH_MAP_G = depth_map_arg
    NBOOT_G = nboot_arg


def process_island(pid):

    sub = ISLAND_TABLES.get(pid)

    if sub is None or sub.empty:
        return None

    D = DEPTH_MAP_G[pid]
    samples = sub["sample"].unique()

    sample_arrays = {}
    sample_pops = {}
    for sample in samples:
        ss = sub[sub["sample"] == sample]
        expanded = np.repeat(ss["pattern"].values, ss["count"].values)
        if len(expanded) >= D:
            sample_arrays[sample] = expanded
            sample_pops[sample] = ss["population"].iloc[0]

    if not sample_arrays:
        return None

    north_vals = []
    south_vals = []
    delta_vals = []

    for b in range(NBOOT_G):

        n_entropy = []
        s_entropy = []

        for sample, expanded in sample_arrays.items():

            sampled = np.random.choice(expanded, size=D, replace=False)
            h = entropy(sampled)
            pop = sample_pops[sample]

            if pop == "north":
                n_entropy.append(h)
            elif pop == "south":
                s_entropy.append(h)

        if len(n_entropy) > 0 and len(s_entropy) > 0:
            n_mean = np.mean(n_entropy)
            s_mean = np.mean(s_entropy)
            north_vals.append(n_mean)
            south_vals.append(s_mean)
            delta_vals.append(n_mean - s_mean)

    if len(delta_vals) == 0:
        return None

    delta = np.mean(delta_vals)
    delta_ci_low = np.percentile(delta_vals, 2.5)
    delta_ci_high = np.percentile(delta_vals, 97.5)
    delta_array = np.array(delta_vals)

    if delta > 0:
        n_extreme = np.sum(delta_array <= 0)
    elif delta < 0:
        n_extreme = np.sum(delta_array >= 0)
    else:
        n_extreme = len(delta_array)

    p_value = (n_extreme + 1) / (len(delta_array) + 1)

    return {
        "cpg_island_id": pid,
        "north_mean": np.mean(north_vals),
        "south_mean": np.mean(south_vals),
        "delta": delta,
        "delta_ci_low": delta_ci_low,
        "delta_ci_high": delta_ci_high,
        "p_value": p_value,
        "north_ci_low": np.percentile(north_vals, 2.5),
        "north_ci_high": np.percentile(north_vals, 97.5),
        "south_ci_low": np.percentile(south_vals, 2.5),
        "south_ci_high": np.percentile(south_vals, 97.5),
    }


print("[STEP] starting parallel bootstrap", flush=True)

all_pids = list(depth_map.keys())
results = []

with Pool(processes=workers, initializer=init_worker,
          initargs=(ISLAND_TABLES, depth_map, n_boot)) as pool:
    for i, r in enumerate(pool.imap_unordered(process_island, all_pids, chunksize=25)):
        if i % 500 == 0:
            print(f"[PROGRESS] {i}/{len(all_pids)}", flush=True)
        if r is not None:
            results.append(r)


# Create output table

out = pd.DataFrame(
    results
)


# FDR correction


p = out["p_value"].to_numpy()

order = np.argsort(p)

p_sorted = p[order]

m = len(p_sorted)


fdr_sorted = (
    p_sorted * m
    / np.arange(
        1,
        m + 1
    )
)


# Enforce monotonicity
fdr_sorted = np.minimum.accumulate(
    fdr_sorted[::-1]
)[::-1]


# FDR cannot exceed 1
fdr_sorted = np.minimum(
    fdr_sorted,
    1.0
)


# Put FDR values back in original row order
out["FDR"] = np.empty(m)

out.loc[
    order,
    "FDR"
] = fdr_sorted


out = out.sort_values(
    "cpg_island_id"
)


# Write output file

out.to_csv(
    out_file,
    sep="\t",
    index=False
)
