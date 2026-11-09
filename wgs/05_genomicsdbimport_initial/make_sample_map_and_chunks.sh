#!/bin/bash
for f in *.g.vcf.gz; do
    s=$(basename "$f" .g.vcf.gz)
    echo -e "${s}\t$(realpath "$f")"
done > cohort.sample_map

bcftools view -h <ONE_SAMPLE>.g.vcf.gz | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/' > contigs.txt

mkdir -p contig_chunks
split -l 100 contigs.txt contig_chunks/chunk_ --additional-suffix=.list

cd contig_chunks
i=1
for f in chunk_*.list; do
  mv "$f" chunk_${i}.list
  i=$((i+1))
done
