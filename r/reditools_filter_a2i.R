#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[[hit + 1]]
}

raw <- get_arg("--raw")
out <- get_arg("--out")
sample <- get_arg("--sample", "sample")
condition <- get_arg("--condition", NA_character_)
strandedness <- tolower(get_arg("--strandedness", "auto"))
min_meanq <- as.numeric(get_arg("--min-meanq", "25"))
min_coverage <- as.numeric(get_arg("--min-coverage", "12"))
min_frequency <- as.numeric(get_arg("--min-frequency", "0.03"))
summary_out <- get_arg("--summary-out", NA_character_)

if (is.null(raw) || is.null(out)) stop("Usage: reditools_filter_a2i.R --raw RAW --out OUT [--sample ID] [--condition COND]")
if (!file.exists(raw)) stop("Missing REDItools raw file: ", raw)

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

dt <- fread(raw, showProgress = FALSE)
required_any <- c("Region", "Position", "AllSubs", "Strand", "Frequency")
missing <- setdiff(required_any, names(dt))
if (length(missing)) stop("REDItools file missing required columns: ", paste(missing, collapse = ", "))

# Different REDItools versions/scripts may use Coverage-q30 or Coverage-q25.
coverage_col <- if ("Coverage-q30" %in% names(dt)) "Coverage-q30" else if ("Coverage-q25" %in% names(dt)) "Coverage-q25" else NA_character_
if (is.na(coverage_col)) stop("REDItools file missing Coverage-q30 or Coverage-q25 column")
if (!("MeanQ" %in% names(dt))) dt[, MeanQ := Inf]

dt <- dt[MeanQ >= min_meanq & get(coverage_col) >= min_coverage & Frequency >= min_frequency]

# A-to-I appears as A>G on the transcript/reference-orientation side and T>C on
# the opposite strand. User's existing scripts used opposite Strand filters for
# reverse-stranded versus non-reverse libraries; make that a parameter here.
if (strandedness %in% c("reverse", "fr-firststrand", "firststrand")) {
  dt <- dt[(AllSubs == "AG" & Strand == 1) | (AllSubs == "TC" & Strand == 2)]
} else if (strandedness %in% c("forward", "fr-secondstrand", "secondstrand", "unstranded", "auto", "fr-unstranded")) {
  dt <- dt[(AllSubs == "AG" & Strand == 2) | (AllSubs == "TC" & Strand == 1)]
}

if (nrow(dt)) setorder(dt, -Frequency)
fwrite(dt, file = out, sep = "\t")

summary_dt <- data.table(
  sample = sample,
  condition = condition,
  n_filtered_sites = nrow(dt),
  mean_frequency = if (nrow(dt)) mean(dt$Frequency, na.rm = TRUE) else NA_real_,
  median_frequency = if (nrow(dt)) median(dt$Frequency, na.rm = TRUE) else NA_real_
)

if (!is.na(summary_out)) {
  dir.create(dirname(summary_out), recursive = TRUE, showWarnings = FALSE)
  fwrite(summary_dt, file = summary_out, sep = "\t")
}
