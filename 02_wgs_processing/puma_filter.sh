#!/bin/bash

# =========================================
# Filter PUMA_allsites.vcf.gz for CpG analysis
# Interactive run (no SLURM)
# =========================================

# Input VCF
VCF="PUMA_allsites.vcf.gz"
OUT_DIR="./vcf_filtered_ps"
mkdir -p "$OUT_DIR"

# Filtering parameters
MAX_CALLED=0.8       # Fraction of individuals called (1 - MAX_MISSING)
MIN_QUAL=30          # Minimum site quality
MIN_DEPTH=5          # Minimum per-sample depth
MAX_DEPTH=30         # Maximum per-sample depth
KEEP_FILE=""         # Optional: file listing individuals to keep

# Output basename
BASENAME=$(basename "$VCF" .vcf.gz)

echo "Processing $VCF on $(hostname)"

# Build vcftools command
CMD="vcftools --gzvcf $VCF --remove-indels --remove-filtered-all \
    --min-alleles 2 --max-alleles 2 \
    --max-missing $MAX_CALLED \
    --minQ $MIN_QUAL \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH \
    --recode --recode-INFO-all --stdout"

# Add --keep if provided
if [[ -f "$KEEP_FILE" ]]; then
    CMD="$CMD --keep $KEEP_FILE"
fi

# Run filtering and pipe through bcftools to remove multiallelics, compress
$CMD | bcftools view -e 'N_ALT>1' -O z -o "$OUT_DIR/${BASENAME}_cpg.vcf.gz"

# Index output
bcftools index "$OUT_DIR/${BASENAME}_cpg.vcf.gz"

# Report number of passing sites
PASS=$(bcftools view -H "$OUT_DIR/${BASENAME}_cpg.vcf.gz" | wc -l)
echo "Sites passing filters for $BASENAME: $PASS"

echo "Filtered VCF written to $OUT_DIR/${BASENAME}_cpg.vcf.gz"

