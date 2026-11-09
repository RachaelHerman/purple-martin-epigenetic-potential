#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(epialleleR)
})

# Set up arguments for running script

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript extract_incremental_v2.R <bam> <bed> <chunk_id> <outdir>")
}

bam <- normalizePath(args[1], mustWork = TRUE)
bed_file <- normalizePath(args[2], mustWork = TRUE)
chunk_id <- args[3]
outdir<- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sample <- tools::file_path_sans_ext(basename(bam))


# Load CpG islands

bed <- fread(bed_file, header = FALSE)
if (ncol(bed) < 3) stop("BED must have at least 3 columns")
setnames(bed, 1:3, c("chrom", "start", "end"))

if (ncol(bed) >= 4) {
  bed[, cpg_island_id := as.character(get(names(bed)[4]))]
} else {
  bed[, cpg_island_id := paste0("cpg_island_", seq_len(.N))]
}
setorder(bed, chrom, start, end)


# Output files (headers)

long_file   <- file.path(outdir, paste0(sample, ".", chunk_id, ".epiallele_long.tsv"))
status_file <- file.path(outdir, paste0(sample, ".", chunk_id, ".epiallele_status.tsv"))

fwrite(data.table(
  cpg_island_id = character(), chrom = character(), status = character(),
  n_patterns = integer(), n_cpg_cols = integer(),
  runtime_sec = numeric(), message = character()
), status_file, sep = "\t")

fwrite(data.table(
  seqnames = character(), strand = character(), start = integer(), end = integer(),
  nbase = integer(), beta = numeric(), pattern = character(),
  position = character(), state = character(),
  cpg_island_id = character(), chrom = character(),
  start_region = integer(), end_region = integer(), sample = character()
), long_file, sep = "\t")

write_status <- function(id, chr, status, n_patterns, n_cpg_cols, t0, msg) {
  fwrite(data.table(
    cpg_island_id = id, chrom = chr, status = status,
    n_patterns = n_patterns, n_cpg_cols = n_cpg_cols,
    runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    message = msg
  ), status_file, sep = "\t", append = TRUE, col.names = FALSE)
}

meta_names <- c("seqnames", "strand", "start", "end", "nbase", "beta", "pattern")


# Extract methylation patterns one CpG island at a time

for (i in seq_len(nrow(bed))) {

  region <- bed[i]
  id  <- region$cpg_island_id
  chr <- region$chrom
  t0  <- Sys.time()

  # Temporary single-island BED
  tmp_bed <- file.path(outdir, paste0("tmp_", chunk_id, "_", i, ".bed"))
  fwrite(region[, .(chrom, start, end)], tmp_bed, sep = "\t", col.names = FALSE)

  ep <- tryCatch(
    extractPatterns(bam = bam, bed = tmp_bed, verbose = FALSE),
    error = function(e) list(.error = e$message)
  )
  invisible(file.remove(tmp_bed))

  if (is.list(ep) && !is.null(ep$.error)) {
    write_status(id, chr, "ERROR", 0, 0, t0, ep$.error)
    next
  }
  if (is.null(ep) || nrow(ep) == 0) {
    write_status(id, chr, "EMPTY", 0, 0, t0, "No epiallele patterns returned")
    next
  }

  dt <- as.data.table(ep)
  meta_cols <- intersect(meta_names, colnames(dt))
  cpg_cols  <- setdiff(colnames(dt), meta_cols)

  if (length(cpg_cols) == 0) {
    write_status(id, chr, "NO_CPG", nrow(dt), 0, t0, "No CpG columns in epiallele output")
    next
  }

  # Long format: one row per pattern x CpG position, empty states dropped
  long_dt <- melt(dt, id.vars = meta_cols, measure.vars = cpg_cols,
                  variable.name = "position", value.name = "state")
  long_dt <- long_dt[!is.na(state) & state != ""]

  if (nrow(long_dt) == 0) {
    write_status(id, chr, "FILTERED_EMPTY", nrow(dt), length(cpg_cols), t0,
                 "All CpG states filtered out")
    next
  }

  long_dt[, `:=`(cpg_island_id = id, chrom = chr,
                 start_region = region$start, end_region = region$end,
                 sample = sample)]
  fwrite(long_dt, long_file, sep = "\t", append = TRUE, col.names = FALSE)

  write_status(id, chr, "PASS", nrow(dt), length(cpg_cols), t0, "OK")
}
