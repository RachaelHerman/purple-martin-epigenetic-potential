[README.md](https://github.com/user-attachments/files/34583438/README.md)
# Whole-genome resequencing pipeline for purple martins (*Progne subis*)

Scripts used to process whole-genome resequencing (WGS) data and call variants for 20 purple martins sampled at the northern and southern edges of the breeding range.

**Author:** Rachael Herman

## Reference genome

*Progne subis* UM_PSub_1 (NCBI GCA_022316685.1)

## Data

Raw reads are available in the NCBI Sequence Read Archive under BioProject [PRJNA1373576](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1373576).

## Software

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

## Repository structure

```
WGS_pipeline.Rmd                 Annotated version of the full pipeline
wgs/
├── 01_split_reads/              Split raw FASTQ files into 1,000,000-read chunks
├── 02_adapter_removal_mapping/  AdapterRemoval and BWA-MEM mapping
├── 03_merge_bams/               Merge BAMs, add read groups, mark duplicates
├── 04_haplotypecaller/          Parallel GATK HaplotypeCaller (GVCF mode)
├── 05_genomicsdbimport_initial/ Sample map, contig chunks, GenomicsDBImport
├── 06_genotypegvcfs_initial/    Joint genotyping (initial round)
├── 07_hard_filter_known_sites/  Hard filtering to make the known-sites set
├── 08_bqsr/                     Base quality score recalibration
├── 10_genomicsdbimport_final/   GenomicsDBImport on recalibrated GVCFs
├── 11_genotypegvcfs_final/      Joint genotyping (final round)
├── 12_merge_final_vcf/          Merge genotyped chunks into one VCF
└── 13_allsites_genotyping/      Joint genotyping including invariant sites
```

## Pipeline steps

| Step | Script(s) | Description |
|---|---|---|
| 1 | `puma_split.sb`, `split_fastq_MP_2.py`, `split_fastq.py` | Split raw reads into chunks for parallel processing |
| 2 | `submit_map.sh`, `merge_map_clean_allsets.py` | Trim adapters and merge overlapping read pairs (AdapterRemoval); map with BWA-MEM |
| 3 | `submit_merge_bams.sh`, `merge_bams_script.sh` | Merge BAMs per sample, add read groups, mark duplicates (Picard) |
| 4 | `hapcaller_batch.sb`, `parallel_hapcaller.py` | Call variants per sample with HaplotypeCaller in GVCF mode |
| 5 | `make_sample_map_and_chunks.sh`, `genomicsDBImport_init.sb` | Build the sample map, split contigs into 30 chunks, import GVCFs |
| 6 | `genotype_chunks_init.sb` | Joint genotyping per contig chunk |
| 7 | `hard_filter_known_sites.sh` | Hard-filter SNPs and indels to make known sites for BQSR |
| 8 | `bqsr.sb` | BaseRecalibrator and ApplyBQSR |
| 9 | Step 4 scripts | HaplotypeCaller rerun on the recalibrated BAMs, with output to `final_var_call/` |
| 10 | `cohort.sample_map`, `genomicDBImport.sb` | GenomicsDBImport of the recalibrated GVCFs |
| 11 | `genotype_chunks.sb` | Final joint genotyping per contig chunk |
| 12 | `merge_final_chunks.sh` | Merge genotyped chunks into the final VCF |
| 13 | `genotype_chunks_allsites.sb` | Joint genotyping with `--include-non-variant-sites` |

## Hard filters (Step 7)

- SNPs removed if QUAL < 50, QD < 2, FS > 60, SOR > 3 or MQ < 40
- Indels removed if QD < 2, QUAL < 50 or FS > 60

## Running the scripts

- The scripts were run on a SLURM cluster. Update the `/path/to/...` placeholders, module names and SLURM resource settings for your system.
- Scripts call helper scripts by file name (for example, `hapcaller_batch.sb` runs `parallel_hapcaller.py`), so run each step from inside its folder, or edit the paths.
- `hapcaller_batch.sb` expects the BAM path as an environment variable: `sbatch --export=BAM=/path/to/sample.bam hapcaller_batch.sb`.

## Note on samples

PM141, PM150 and PM223 (*P. subis*, a different subspecies) were included in variant calling and joint genotyping for a separate PCA. They were removed before variant filtering and were not used in any analyses reported in the manuscript.

## License

MIT License. See [LICENSE](LICENSE).
