#!/bin/bash
#SBATCH --job-name=PM141_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --mem=1000G
#SBATCH --time=48:00:00
#SBATCH --partition=delmore_lab1
#SBATCH --account=delmore_lab
#SBATCH --output=$LOG_DIR/%A_%a.merge.out
#SBATCH --error=$LOG_DIR/%A_%a.merge.err


# Usage: sbatch merge_align_pipeline.sh SAMPLE_NAME
SAMPLE=$1
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    exit 1
fi

# Directories
BAMS_DIR="/insomnia001/depts/delmore_lab/users/rwh2140/puma_wgs/${SAMPLE}/${SAMPLE}/adaptrem/bams"
FASTQ_DIR="/insomnia001/depts/delmore_lab/users/rwh2140/puma_wgs/${SAMPLE}/${SAMPLE}/adaptrem/fastq"
OUT_DIR="$BAMS_DIR/final"
mkdir -p "$OUT_DIR"
cd "$BAMS_DIR" || { echo "BAM dir not found: $BAMS_DIR"; exit 1; }

# Tools
BWA="/insomnia001/depts/delmore_lab/projects/tools/bwa-mem2/bwa-mem2.avx2"
SAMTOOLS="samtools"
JAVA="/insomnia001/depts/delmore_lab/projects/tools/openjdk17/bin/java"
PICARD="/insomnia001/depts/delmore_lab/projects/tools/picard/picard.jar"

THREADS=160

# Function to merge split BAMs (if they exist) and add read groups
process_bam_type() {
    local TYPE=$1
    SPLIT_BAMS=( $(ls ${SAMPLE}_*_split_*_${TYPE}.bam 2>/dev/null) )
    if [ ${#SPLIT_BAMS[@]} -eq 0 ]; then
        echo "No $TYPE BAMs found, skipping."
        return
    fi

    MERGED_BAM="$OUT_DIR/${SAMPLE}_${TYPE}_merged.bam"
    echo "Merging ${#SPLIT_BAMS[@]} $TYPE BAMs..."
    $SAMTOOLS merge -@ $THREADS "$MERGED_BAM" "${SPLIT_BAMS[@]}"

    SORTED_BAM="$OUT_DIR/${SAMPLE}_${TYPE}_merged.sort.bam"
    echo "Sorting $MERGED_BAM..."
    $SAMTOOLS sort -@ $THREADS -o "$SORTED_BAM" "$MERGED_BAM"
    $SAMTOOLS index "$SORTED_BAM"

    RG_BAM="$OUT_DIR/${SAMPLE}_${TYPE}_merged.sort.RG.bam"
    echo "Adding read groups..."
    $JAVA -Xmx500G -jar $PICARD AddOrReplaceReadGroups \
        I="$SORTED_BAM" \
        O="$RG_BAM" \
        LB=Veeramah PL=COMPLETE PU=BGI SM="$SAMPLE" CN=BGI RGID="${SAMPLE}_${TYPE}" \
        VALIDATION_STRINGENCY=SILENT

    FINAL_BAM="$OUT_DIR/${SAMPLE}_${TYPE}_final.Mkdup.bam"
    echo "Marking duplicates..."
    $JAVA -Xmx500G -jar $PICARD MarkDuplicates \
        I="$RG_BAM" \
        O="$FINAL_BAM" \
        M="$OUT_DIR/${SAMPLE}_${TYPE}.metrics" \
        AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT
    $SAMTOOLS index "$FINAL_BAM"

    $SAMTOOLS flagstat "$FINAL_BAM" > "${FINAL_BAM}.flagstat"
    echo "$TYPE BAM processing complete: $FINAL_BAM"

    # Clean intermediate files
    rm -f "$MERGED_BAM" "$SORTED_BAM" "$RG_BAM"
}

# Process PE and ME BAMs
process_bam_type "PE"
process_bam_type "ME"

# Merge PE+ME final BAMs if both exist
PE_FINAL="$OUT_DIR/${SAMPLE}_PE_final.Mkdup.bam"
ME_FINAL="$OUT_DIR/${SAMPLE}_ME_final.Mkdup.bam"
if [ -f "$PE_FINAL" ] && [ -f "$ME_FINAL" ]; then
    FINAL_BAM="$OUT_DIR/${SAMPLE}_PE_ME_final.Mkdup.bam"
    echo "Merging PE + ME BAMs..."
    $SAMTOOLS merge -@ $THREADS "$FINAL_BAM" "$PE_FINAL" "$ME_FINAL"

    SORTED_FINAL="$OUT_DIR/${SAMPLE}_PE_ME_final.sort.bam"
    $SAMTOOLS sort -@ $THREADS -o "$SORTED_FINAL" "$FINAL_BAM"
    $SAMTOOLS index "$SORTED_FINAL"

    $JAVA -Xmx500G -jar $PICARD MarkDuplicates \
        I="$SORTED_FINAL" \
        O="$FINAL_BAM" \
        M="$OUT_DIR/${SAMPLE}_PE_ME.metrics" \
        AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT
    $SAMTOOLS index "$FINAL_BAM"

    $SAMTOOLS flagstat "$FINAL_BAM" > "${FINAL_BAM}.flagstat"
    echo "Final PE+ME BAM complete: $FINAL_BAM"

    # Optional cleanup
    rm -f "$SORTED_FINAL"
fi

echo "All done for sample $SAMPLE!"
