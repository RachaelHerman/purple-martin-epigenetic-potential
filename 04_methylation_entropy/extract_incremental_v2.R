#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(epialleleR)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript extract_incremental.R <bam> <bed> <chunk_id> <outdir>")
}

bam <- normalizePath(args[1], mustWork = TRUE)
bed_file <- normalizePath(args[2], mustWork = TRUE)
chunk_id <- args[3]
outdir <- args[4]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sample <- tools::file_path_sans_ext(basename(bam))

message("=== START ===")
message("Sample: ", sample)
message("Chunk: ", chunk_id)

# ----------------------------
# Load BED
# ----------------------------
bed <- fread(bed_file, header = FALSE)

if (ncol(bed) < 3) stop("BED must have at least 3 columns")

setnames(bed, 1:3, c("chrom", "start", "end"))

if (ncol(bed) >= 4) {
  bed[, promoter_id := as.character(get(names(bed)[4]))]
} else {
  bed[, promoter_id := paste0("promoter_", seq_len(.N))]
}

setorder(bed, chrom, start, end)

# ----------------------------
# OUTPUT FILES
# ----------------------------
long_file <- file.path(outdir, paste0(sample, ".", chunk_id, ".epiallele_long.tsv"))
status_file <- file.path(outdir, paste0(sample, ".", chunk_id, ".epiallele_status.tsv"))

# status header
fwrite(data.table(
  promoter_id=character(),
  chrom=character(),
  status=character(),
  n_patterns=integer(),
  n_cpg_cols=integer(),
  runtime_sec=numeric(),
  message=character()
), status_file, sep="\t")

# long header
fwrite(data.table(
  seqnames=character(),
  strand=character(),
  start=integer(),
  end=integer(),
  nbase=integer(),
  beta=numeric(),
  pattern=character(),
  position=character(),
  state=character(),
  promoter_id=character(),
  chrom=character(),
  start_region=integer(),
  end_region=integer(),
  sample=character()
), long_file, sep="\t")

# ----------------------------
# MAIN LOOP
# ----------------------------
message("Processing promoters: ", nrow(bed))
flush.console()

for (i in seq_len(nrow(bed))) {

  region <- bed[i]
  pid <- region$promoter_id
  chr <- region$chrom

  message("START ", i, " ", pid)
  flush.console()

  t0 <- Sys.time()

  # ----------------------------
  # UNIQUE TEMP BED PER PROMOTER
  # ----------------------------
  tmp_bed <- file.path(
    outdir,
    paste0("tmp_", chunk_id, "_", i, ".bed")
  )

  fwrite(region[, .(chrom, start, end)], tmp_bed,
         sep = "\t", col.names = FALSE)

  ep <- tryCatch(
    extractPatterns(bam = bam, bed = tmp_bed, verbose = FALSE),
    error = function(e) list(.error = e$message)
  )

  # optional cleanup (safe)
  file.remove(tmp_bed)

  # ----------------------------
  # ERROR
  # ----------------------------
  if (is.list(ep) && !is.null(ep$.error)) {

    fwrite(data.table(
      promoter_id=pid,
      chrom=chr,
      status="ERROR",
      n_patterns=0,
      n_cpg_cols=0,
      runtime_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
      message=ep$.error
    ), status_file, sep="\t", append=TRUE, col.names=FALSE)

    message("ERROR ", pid, ": ", ep$.error)
    next
  }

  # ----------------------------
  # EMPTY
  # ----------------------------
  if (is.null(ep) || nrow(ep) == 0) {

    fwrite(data.table(
      promoter_id=pid,
      chrom=chr,
      status="EMPTY",
      n_patterns=0,
      n_cpg_cols=0,
      runtime_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
      message="No epiallele patterns returned"
    ), status_file, sep="\t", append=TRUE, col.names=FALSE)

    message("EMPTY ", pid)
    next
  }

  dt <- as.data.table(ep)

  meta_cols <- intersect(
    c("seqnames","strand","start","end","nbase","beta","pattern"),
    colnames(dt)
  )

  cpg_cols <- setdiff(colnames(dt), meta_cols)

  # ----------------------------
  # NO CpG STRUCTURE
  # ----------------------------
  if (length(cpg_cols) == 0) {

    fwrite(data.table(
      promoter_id=pid,
      chrom=chr,
      status="NO_CPG",
      n_patterns=nrow(dt),
      n_cpg_cols=0,
      runtime_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
      message="No CpG columns in epiallele output"
    ), status_file, sep="\t", append=TRUE, col.names=FALSE)

    message("NO_CPG ", pid)
    next
  }

  # ----------------------------
  # LONG FORMAT
  # ----------------------------
  long_dt <- melt(
    dt,
    id.vars = meta_cols,
    measure.vars = cpg_cols,
    variable.name = "position",
    value.name = "state"
  )

  long_dt <- long_dt[!is.na(state) & state != ""]

  if (nrow(long_dt) == 0) {

    fwrite(data.table(
      promoter_id=pid,
      chrom=chr,
      status="FILTERED_EMPTY",
      n_patterns=nrow(dt),
      n_cpg_cols=length(cpg_cols),
      runtime_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
      message="All CpG states filtered out"
    ), status_file, sep="\t", append=TRUE, col.names=FALSE)

    message("FILTERED_EMPTY ", pid)
    next
  }

  # annotate
  long_dt[, `:=`(
    promoter_id = pid,
    chrom = chr,
    start_region = region$start,
    end_region = region$end,
    sample = sample
  )]

  fwrite(long_dt, long_file, sep="\t", append=TRUE, col.names=FALSE)

  fwrite(data.table(
    promoter_id=pid,
    chrom=chr,
    status="PASS",
    n_patterns=nrow(dt),
    n_cpg_cols=length(cpg_cols),
    runtime_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
    message="OK"
  ), status_file, sep="\t", append=TRUE, col.names=FALSE)

  message("PASS ", pid, " patterns=", nrow(dt))

  flush.console()
}

message("=== DONE ===")
