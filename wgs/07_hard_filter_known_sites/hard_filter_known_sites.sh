#!/bin/bash
ls genotyped_chunks/chunk_*.vcf.gz | sort -V > chunks.list

gatk GatherVcfs -I chunks.list -O merged_genotyped_chunks.vcf.gz
gatk IndexFeatureFile -I merged_genotyped_chunks.vcf.gz

bcftools view -v snps merged_genotyped_chunks.vcf.gz -Oz -o all_snps.vcf.gz
bcftools index -f all_snps.vcf.gz

bcftools view -e 'QUAL<50 || INFO/QD<2 || INFO/FS>60 || INFO/SOR>3 || INFO/MQ<40' \
  all_snps.vcf.gz -Oz -o known_sites.snps.vcf.gz
bcftools index -f known_sites.snps.vcf.gz

bcftools view -v indels merged_genotyped_chunks.vcf.gz -Oz -o all_indels.vcf.gz
bcftools index -f all_indels.vcf.gz

bcftools view -e 'INFO/QD<2 || QUAL<50 || INFO/FS>60' \
  all_indels.vcf.gz -Oz -o known_sites.indels.vcf.gz
bcftools index -f known_sites.indels.vcf.gz

bcftools concat known_sites.snps.vcf.gz known_sites.indels.vcf.gz \
  -a -O z -o known_sites.snp_indel.vcf.gz
