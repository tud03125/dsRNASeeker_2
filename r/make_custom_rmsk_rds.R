#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(GenomeInfoDb)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop(
        "Usage: Rscript make_custom_rmsk_rds.R ",
        "<TE_only.bed> <output.rds>"
    )
}

bed_file <- args[[1]]
output_rds <- args[[2]]

if (!file.exists(bed_file)) {
    stop("Input BED does not exist: ", bed_file)
}

if (file.info(bed_file)$size == 0) {
    stop("Input BED is empty: ", bed_file)
}

message("Reading TE-only BED: ", bed_file)

bed <- read.delim(
    bed_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    col.names = c(
        "seqnames",
        "start0",
        "end",
        "TE_id",
        "score",
        "strand",
        "repName",
        "repClass",
        "repFamily"
    )
)

if (nrow(bed) == 0) {
    stop("No records were read from the TE-only BED.")
}

bed$start0 <- suppressWarnings(as.integer(bed$start0))
bed$end <- suppressWarnings(as.integer(bed$end))

valid <- (
    !is.na(bed$seqnames) &
    !is.na(bed$start0) &
    !is.na(bed$end) &
    bed$start0 >= 0 &
    bed$end > bed$start0 &
    bed$strand %in% c("+", "-") &
    !is.na(bed$TE_id) &
    bed$TE_id != ""
)

if (!all(valid)) {
    message("Removing invalid BED rows: ", sum(!valid))
    bed <- bed[valid, , drop = FALSE]
}

if (nrow(bed) == 0) {
    stop("No valid rows remained after BED validation.")
}

message("Constructing GRanges from ", format(nrow(bed), big.mark=","), " records")

rmsk <- GRanges(
    seqnames = bed$seqnames,
    ranges = IRanges(
        start = bed$start0 + 1L,
        end = bed$end
    ),
    strand = bed$strand
)

unique_ids <- make.unique(as.character(bed$TE_id))

names(rmsk) <- unique_ids

mcols(rmsk)$TE_id <- unique_ids
mcols(rmsk)$repName <- as.character(bed$repName)
mcols(rmsk)$repClass <- as.character(bed$repClass)
mcols(rmsk)$repFamily <- as.character(bed$repFamily)

# Retained for compatibility with downstream annotation code.
mcols(rmsk)$repStart <- NA_integer_
mcols(rmsk)$repEnd <- NA_integer_

dir.create(
    dirname(output_rds),
    recursive = TRUE,
    showWarnings = FALSE
)

message("Saving GRanges RDS: ", output_rds)

saveRDS(
    rmsk,
    file = output_rds,
    compress = TRUE
)

cat("[OK] class:", class(rmsk)[1], "\n")
cat("[OK] ranges:", length(rmsk), "\n")
cat("[OK] seqlevels:", length(seqlevels(rmsk)), "\n")
cat(
    "[OK] metadata:",
    paste(names(mcols(rmsk)), collapse = ","),
    "\n"
)
cat("[OK] saved:", output_rds, "\n")
