#!/usr/bin/env python3
import os
import sys
from subprocess import Popen
import math

# ==== Configure your tools ====
JAVA = "/insomnia001/depts/delmore_lab/projects/tools/openjdk17/bin/java"
PICARD = "/insomnia001/depts/delmore_lab/projects/tools/picard/picard.jar"
SAMTOOLS_THREADS = 160
BATCH_SIZE = 8  # number of BAMs to merge at once

def run(cmd):
    print(f"Running: {cmd}")
    ret = Popen(cmd, shell=True)
    ret.wait()
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")

def list_bams(sample_dir, sample_name, bam_type):
    """List BAM files for a given type (PE or ME), ignoring MEt"""
    bams_dir = os.path.join(sample_dir, "adaptrem", "bams")
    files = [f for f in os.listdir(bams_dir)
             if f.startswith(sample_name) and f.endswith(f"_{bam_type}.bam")]
    # Sort by split number (assumes format: SAMPLE_1_split_XXX_PE.bam)
    files.sort(key=lambda x: int(x.split('_')[3]))
    return [os.path.join(bams_dir, f) for f in files]

def merge_batch(bam_list, prefix, batch_num):
    """Merge a batch of BAMs and return the merged BAM path"""
    merged = f"{prefix}_batch{batch_num}.bam"
    run(f"samtools merge -@ {SAMTOOLS_THREADS} {merged} {' '.join(bam_list)}")
    sorted_bam = merged.replace(".bam", ".sort.bam")
    run(f"samtools sort -@ {SAMTOOLS_THREADS} -o {sorted_bam} {merged}")
    run(f"samtools index {sorted_bam}")
    return sorted_bam

def merge_and_mark_all(bam_list, out_prefix):
    """Merge BAMs in batches, then final merge + mark duplicates"""
    if len(bam_list) == 0:
        raise ValueError("No BAMs found to merge!")

    # Step 1: Batch merge
    batch_files = []
    total_batches = math.ceil(len(bam_list) / BATCH_SIZE)
    for i in range(total_batches):
        batch = bam_list[i*BATCH_SIZE:(i+1)*BATCH_SIZE]
        batch_sorted = merge_batch(batch, out_prefix, i+1)
        batch_files.append(batch_sorted)

    # Step 2: Final merge of batches
    final_merged = f"{out_prefix}.merged.bam"
    run(f"samtools merge -@ {SAMTOOLS_THREADS} {final_merged} {' '.join(batch_files)}")
    final_sorted = final_merged.replace(".bam", ".sort.bam")
    run(f"samtools sort -@ {SAMTOOLS_THREADS} -o {final_sorted} {final_merged}")
    run(f"samtools index {final_sorted}")

    # Step 3: MarkDuplicates
    mkdup_bam = final_sorted.replace(".sort.bam", ".Mkdup.bam")
    metrics = final_sorted.replace(".sort.bam", ".Mkdup.metrics")
    run(f"{JAVA} -jar {PICARD} MarkDuplicates I={final_sorted} O={mkdup_bam} M={metrics} "
        f"REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT AS=TRUE")
    run(f"samtools index {mkdup_bam}")

    return mkdup_bam

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 batch_merge_bams.py SAMPLE_NAME SAMPLE_DIR")
        sys.exit(1)

    sample_name = sys.argv[1]
    sample_dir = sys.argv[2]

    # ==== List PE and ME BAMs ====
    pe_bams = list_bams(sample_dir, sample_name, "PE")
    me_bams = list_bams(sample_dir, sample_name, "ME")

    # ==== Merge and mark duplicates ====
    print("Merging PE BAMs...")
    pe_final = merge_and_mark_all(pe_bams, os.path.join(sample_dir, "adaptrem", "bams", f"{sample_name}_PE"))

    print("Merging ME BAMs...")
    me_final = merge_and_mark_all(me_bams, os.path.join(sample_dir, "adaptrem", "bams", f"{sample_name}_ME"))

    # ==== Merge PE + ME ====
    print("Merging PE + ME BAMs...")
    final_bam = merge_and_mark_all([pe_final, me_final], os.path.join(sample_dir, "adaptrem", "bams", f"{sample_name}_PE_ME"))

    print(f"Final BAM: {final_bam}")
    print("All done!")

if __name__ == "__main__":
    main()

