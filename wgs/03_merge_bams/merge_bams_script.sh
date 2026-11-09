#!/bin/bash
# Usage: ./merge_bams_script.sh SAMPLE_NAME

SAMPLE=$1

BAMS_DIR="/path/to/puma_wgs/${SAMPLE}/${SAMPLE}/adaptrem/bams"
JAVA="/path/to/openjdk17/bin/java"
PICARD="/path/to/picard.jar"

cd "$BAMS_DIR" || exit 1

process_bam_type() {
    local TYPE=$1

    SPLIT_BAMS=( $(ls ${SAMPLE}_*_split_*_${TYPE}.bam 2>/dev/null) )
    if [ ${#SPLIT_BAMS[@]} -eq 0 ]; then
        return
    fi

    MERGED_BAM="${SAMPLE}_${TYPE}_merged.bam"
    samtools merge "$MERGED_BAM" "${SPLIT_BAMS[@]}"

    SORTED_BAM="${SAMPLE}_${TYPE}_merged.sort.bam"
    samtools sort -o "$SORTED_BAM" "$MERGED_BAM"
    samtools index "$SORTED_BAM"

    RG_BAM="${SAMPLE}_${TYPE}_merged.sort.RG.bam"
    $JAVA -jar $PICARD AddOrReplaceReadGroups \
        I="$SORTED_BAM" \
        O="$RG_BAM" \
        LB=Veeramah PL=COMPLETE PU=BGI SM="$SAMPLE" CN=BGI RGID="${SAMPLE}_${TYPE}" \
        VALIDATION_STRINGENCY=SILENT

    FINAL_SORT_BAM="${SAMPLE}_${TYPE}_merged.sort.RG.final.bam"
    samtools sort -o "$FINAL_SORT_BAM" "$RG_BAM"
    samtools index "$FINAL_SORT_BAM"

    MKDUP_BAM="${SAMPLE}_${TYPE}_merged.sort.RG.final.Mkdup.bam"
    $JAVA -jar $PICARD MarkDuplicates \
        I="$FINAL_SORT_BAM" \
        O="$MKDUP_BAM" \
        M="${SAMPLE}_${TYPE}_merged.metrics" \
        AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT
    samtools index "$MKDUP_BAM"

    samtools flagstat "$MKDUP_BAM" > "${MKDUP_BAM}.flagstat"

    rm -f "$MERGED_BAM" "$SORTED_BAM" "$RG_BAM" "$FINAL_SORT_BAM"
}

process_bam_type "PE"
process_bam_type "ME"

PE_FINAL="${SAMPLE}_PE_merged.sort.RG.final.Mkdup.bam"
ME_FINAL="${SAMPLE}_ME_merged.sort.RG.final.Mkdup.bam"

if [ -f "$PE_FINAL" ] && [ -f "$ME_FINAL" ]; then
    FINAL_BAM="${SAMPLE}_PE_ME_merged.sort.Mkdup.bam"
    samtools merge "$FINAL_BAM" "$PE_FINAL" "$ME_FINAL"

    SORTED_FINAL="${FINAL_BAM%.bam}.sort.bam"
    samtools sort -o "$SORTED_FINAL" "$FINAL_BAM"
    samtools index "$SORTED_FINAL"

    $JAVA -jar $PICARD MarkDuplicates \
        I="$SORTED_FINAL" \
        O="$FINAL_BAM" \
        M="${SAMPLE}_PE_ME_merged.metrics" \
        AS=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT
    samtools index "$FINAL_BAM"

    samtools flagstat "$FINAL_BAM" > "${FINAL_BAM}.flagstat"
fi
