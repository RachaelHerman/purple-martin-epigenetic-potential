# Genomic and methylation pipelines for purple martins (*Progne subis*)

Scripts used to process whole-genome resequencing (WGS) and enzymatic methyl-seq (EM-seq) data, and to estimate realized methylation entropy of CpG islands, for 20 purple martins sampled at the northern and southern edges of the breeding range.

**Author:** Rachael Herman

## Reference genome

*Progne subis* UM_PSub_1 (NCBI GCA_022316685.1)

## Data

Raw reads are available in the NCBI Sequence Read Archive under BioProject [PRJNA1373576](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1373576).

## Repository structure

```
WGS_pipeline.Rmd                 Annotated WGS pipeline
Methylation_pipeline.Rmd         Annotated EM-seq and entropy pipeline
wgs/                             WGS read processing and variant calling
methylation/                     EM-seq processing and methylation entropy
```

Each pipeline folder contains numbered step folders with the scripts for that step. Each `.Rmd` file describes the same steps with the code shown inline.

---

## WGS pipeline (`wgs/`)

### Software

| Tool | Version |
|---|---|
| AdapterRemoval | 2.3.4 |
| BWA | 0.7.17  |
| Samtools | 1.7 |
| Picard | 3.4.0 |
| GATK | 4.6.2.0 |
| OpenJDK | 17 |
| bcftools | 1.21 |
| Python | 3.13 |
| Biopython | 1.85 |

### Steps

| Step | Folder | Script(s) | Description |
|---|---|---|---|
| 1 | `01_split_reads` | `puma_split.sb`, `split_fastq_MP_2.py`, `split_fastq.py` | Split raw reads into 1,000,000-read chunks for parallel processing |
| 2 | `02_adapter_removal_mapping` | `submit_map.sh`, `merge_map_clean_allsets.py` | Trim adapters and merge overlapping read pairs (AdapterRemoval); map with BWA-MEM |
| 3 | `03_merge_bams` | `submit_merge_bams.sh`, `merge_bams_script.sh` | Merge BAMs per sample, add read groups, mark duplicates (Picard) |
| 4 | `04_haplotypecaller` | `hapcaller_batch.sb`, `parallel_hapcaller.py` | Per-sample variant calling with HaplotypeCaller in GVCF mode |
| 5 | `05_genomicsdbimport_initial` | `make_sample_map_and_chunks.sh`, `genomicsDBImport_init.sb` | Build the sample map, split contigs into 30 chunks, import GVCFs |
| 6 | `06_genotypegvcfs_initial` | `genotype_chunks_init.sb` | Joint genotyping per contig chunk |
| 7 | `07_hard_filter_known_sites` | `hard_filter_known_sites.sh` | Hard-filter SNPs and indels to make known sites for BQSR |
| 8 | `08_bqsr` | `bqsr.sb` | BaseRecalibrator and ApplyBQSR |
| 9 | (Step 4 scripts) | `hapcaller_batch.sb`, `parallel_hapcaller.py` | HaplotypeCaller rerun on the recalibrated BAMs, with output to `final_var_call/` |
| 10 | `10_genomicsdbimport_final` | `cohort.sample_map`, `genomicDBImport.sb` | GenomicsDBImport of the recalibrated GVCFs |
| 11 | `11_genotypegvcfs_final` | `genotype_chunks.sb` | Final joint genotyping per contig chunk |
| 12 | `12_merge_final_vcf` | `merge_final_chunks.sh` | Merge genotyped chunks into the final VCF |
| 13 | `13_allsites_genotyping` | `genotype_chunks_allsites.sb` | Joint genotyping with `--include-non-variant-sites` |

Hard filters (Step 7):

- SNPs removed if QUAL < 50, QD < 2, FS > 60, SOR > 3 or MQ < 40
- Indels removed if QD < 2, QUAL < 50 or FS > 60

PM141, PM150 and PM223 (*P. subis*, a different subspecies) were included in variant calling and joint genotyping for a separate PCA. They were removed before variant filtering and were not used in any analyses reported in the manuscript.

---

## Methylation pipeline (`methylation/`)

### Software

| Tool | Version |
|---|---|
| Bismark | 0.25.1 |
| Bowtie 2 | 2.5.4 |
| Trim Galore | 0.6.10 |
| Cutadapt | 5.2 |
| FastQC | 0.12.1 |
| Samtools | 1.17 (alignment pipeline), 1.22 (CpG island subsetting) |
| Picard | 3.4.0 |
| EMBOSS | [VERSION] |
| R | 4.3.3 (in container) |
| Bioconductor | 3.18 (in container) |
| epialleleR | [VERSION] |
| Python / pandas / NumPy | [VERSION] |
| Apptainer | [VERSION] |

### Steps

| Step | Folder | Script(s) | Description |
|---|---|---|---|
| 1 | `01_bismark_genome_prep` | `bismark_genome_prep.sh` | In-silico bisulfite conversion and Bowtie 2 indexing of the reference |
| 2 | `02_trim_reads` | `trim_emseq.sb` | Adapter and quality trimming (Trim Galore), FastQC on trimmed reads |
| 3 | `03_align_dedup` | `bismark_align_dedup.sb` | Bismark alignment, sorting, read groups, duplicate removal, name-sorting |
| 4 | `04_cpg_islands` | `newcpgreport.sh` | Identify CpG islands with EMBOSS newcpgreport (default settings) |
| 5 | `05_subset_bams_to_cgis` | `subset_bams_to_cgis.sh` | Extract reads overlapping CpG islands, fix mate information, name-sort |
| 6 | `06_epialleler_container` | `epialleleR.def`, `build_container.sh` | Apptainer container for epialleleR |
| 7 | `07_extract_patterns` | `epialleleR_batch_chunk.sb`, `extract_incremental_v2.R` | Extract read-level methylation patterns per CpG island |
| 8 | `08_count_patterns` | `beta_counter.sb`, `beta_epiallele_cgi.py` | Count reads per methylation pattern per sample and CpG island |
| 9 | `09_entropy_subsampling` | `bootstrap_cgi.sb`, `epi_entropy_bootstrap_cgi_v2.py` | Realized methylation entropy (Shannon) by repeated subsampling to a fixed depth per island; North vs South |

---

## Running the scripts

- The scripts were run on a SLURM cluster. Update the `/path/to/...` placeholders, module and conda environment names, and SLURM resource settings for your system.
- Scripts call helper scripts by file name, so run each step from inside its folder or edit the paths.
- `wgs/04_haplotypecaller/hapcaller_batch.sb` expects the BAM path as an environment variable: `sbatch --export=BAM=/path/to/sample.bam hapcaller_batch.sb`.
- `methylation/05_subset_bams_to_cgis/subset_bams_to_cgis.sh` shows the commands for one sample; set `sample` first or loop over samples.
- For `methylation/07_extract_patterns`, the CpG island BED file was split into chunks for parallel processing. Chunking does not affect results, and the full BED file can be passed instead.

## License

MIT License. See [LICENSE](LICENSE).

