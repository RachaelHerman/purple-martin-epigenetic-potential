#!/bin/bash
# Coordinate-sort
samtools sort -o ${sample}_PE_merged.sort.readname.RG.Mkdup.bismark.coord.bam \
    ${sample}_PE_merged.sort.readname.RG.Mkdup.bismark.bam

# Extract reads overlapping CpG islands
samtools view -b -L PUMA_CpG_islands.sorted.bed \
    ${sample}_PE_merged.sort.readname.RG.Mkdup.bismark.coord.bam \
    > ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.coord.bam

# Name-sort, fix mate information, and name-sort again
samtools sort -n -o ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.name.bam \
    ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.coord.bam

samtools fixmate -m ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.name.bam \
    ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.fixmate.bam

samtools sort -n -o ${sample}_PE_merged.sort.readname.RG.Mkdup.CpG_islands.name_sorted.bam \
    ${sample}_PE_merged.sort.readname.RG.Mkdup.islands.fixmate.bam
