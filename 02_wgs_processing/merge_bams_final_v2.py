#!/usr/bin/env python3
import os
import sys
from subprocess import Popen

# Paths to Java and Picard
JAVA = "/insomnia001/depts/delmore_lab/projects/tools/openjdk17/bin/java"
PICARD = "/insomnia001/depts/delmore_lab/projects/tools/picard/picard.jar"

# Number of threads to use for samtools
THREADS = 160

def run(cmd):
    print(f"Running: {cmd}")
    ret = Popen(cmd, shell=True)
    ret.wait()
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")

def list_bams(bam_dir, sample_name, bam_type):
    """List split BAMs for a given type (PE or ME)"""
    files = [f for f in os.listdir(bam_dir)
             if f.startswith(sample_name) and "_split_" in f and f.endswith(f"_{bam_type}.bam")]
    files.sort(key=lambda x: int(x.split('_')[3]))  # sort by split number
    return [os.path.join(bam_dir, f) for f in files]

def merge_sort_rg_mkdup(bam_list, output_prefix):
    """Merge BAMs, sort, add read groups, mark duplicates"""
    if len(bam_list) == 0:
        return None

    merged_bam = f"{output_prefix}.merged.bam"
    sorted_bam = f"{output_prefix}.sort.bam"
    rg_bam = f"{output_prefix}.RG.bam"
    final_bam = f"{output_prefix}.Mkdup.bam"
    metrics_file = f"{output_prefix}.metrics"

    # Merge
    run(f"samtools merge -@ {THREADS} {merged_bam} {' '.join(bam_list)}")

    # Sort
    run(f"samtools sort -@ {THREADS} -o {sorted_bam} {merged_bam}")
    run(f"samtools index {sorted_bam}")

    # Add read groups
    run(f"{JAVA} -jar {PICARD} AddOrReplaceReadGroups I={sorted_bam} O={rg_bam} "
        f"LB=Veeramah PL=COMPLETE PU=BGI SM={output_prefix} CN=BGI RGID={output_prefix} "
        f"VALIDATION_STRINGENCY=SILENT")
    run(f"samtools sort -@ {THREADS} -o {sorted_bam} {rg_bam}")
    run(f"samtools index {sorted_bam}")

    # Mark duplicates
    run(f"{JAVA} -jar {PICARD} MarkDuplicates I={sorted_bam} O={final_bam} M={metrics_file} "
        f"AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT")
    run(f"samtools index {final_bam}")

    # Flagstat
    run(f"samtools flagstat {final_bam} > {final_bam}.flagstat")

    # Cleanup
    for f in [merged_bam, rg_bam, sorted_bam]:
        if os.path.exists(f):
            os.remove(f)

    return final_bam

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 merge_bams_final_fixed.py SAMPLE_NAME BAM_DIR")
        sys.exit(1)

    sample_name = sys.argv[1]
    bam_dir = sys.argv[2]

    if not os.path.exists(bam_dir):
        raise FileNotFoundError(f"BAM directory not found: {bam_dir}")

    # PE BAMs
    pe_bams = list_bams(bam_dir, sample_name, "PE")
    print(f"Found {len(pe_bams)} PE BAMs")
    pe_final = merge_sort_rg_mkdup(pe_bams, os.path.join(bam_dir, f"{sample_name}_PE"))

    # ME BAMs
    me_bams = list_bams(bam_dir, sample_name, "ME")
    print(f"Found {len(me_bams)} ME BAMs")
    me_final = merge_sort_rg_mkdup(me_bams, os.path.join(bam_dir, f"{sample_name}_ME"))

    # Merge PE + ME final BAMs
    if pe_final and me_final:
        final_merged = os.path.join(bam_dir, f"{sample_name}_PE_ME.Mkdup.bam")
        run(f"samtools merge -@ {THREADS} {final_merged} {pe_final} {me_final}")
        run(f"samtools sort -@ {THREADS} -o {final_merged} {final_merged}")
        run(f"{JAVA} -jar {PICARD} MarkDuplicates I={final_merged} O={final_merged} "
            f"M={bam_dir}/{sample_name}_PE_ME.metrics AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT")
        run(f"samtools index {final_merged}")
        run(f"samtools flagstat {final_merged} > {final_merged}.flagstat")
        print(f"Final PE+ME BAM complete: {final_merged}")

if __name__ == "__main__":
    main()

