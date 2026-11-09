#!/bin/bash
ls genotyped_chunks_final/chunk_*.vcf.gz | sort -V > final_chunks.list

gatk GatherVcfs -I final_chunks.list -O final_genotyped.vcf.gz
gatk IndexFeatureFile -I final_genotyped.vcf.gz
