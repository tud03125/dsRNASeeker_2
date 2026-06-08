#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(data.table)
  library(txdbmaker)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx+1]]
}
input_bed <- get_arg("--input-bed", "mouse_te.bed")
gtf_path  <- get_arg("--gtf", NULL)
output_tsv <- get_arg("--output", "mouse_te_txmap.tsv")
if (is.null(gtf_path)) stop("--gtf is required")

te <- fread(input_bed, header = FALSE)
if (ncol(te) < 6) stop("Input BED must have at least 6 columns")
setnames(te, 1:min(10, ncol(te)), c("chrom","start","end","name","score","strand","repFamily","repName","SYMBOL","annotation")[1:min(10, ncol(te))])

gr <- GRanges(seqnames = te$chrom,
              ranges   = IRanges(as.integer(te$start) + 1L, as.integer(te$end)),
              strand   = te$strand)

#txdb <- makeTxDbFromGFF(gtf_path)
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)
tx <- transcripts(txdb)
ov <- findOverlaps(gr, tx, ignore.strand = FALSE)

m <- data.table(
  TE_name   = te$name[queryHits(ov)],
  TE_strand = as.character(strand(gr))[queryHits(ov)],
  TX_strand = as.character(strand(tx))[subjectHits(ov)]
)

m_one <- m[, .(TX_strand = if (uniqueN(TX_strand)==1) TX_strand[1] else NA_character_), by = .(TE_name, TE_strand)]
fwrite(m_one, output_tsv, sep = "	")
message("Wrote: ", output_tsv, " (", nrow(m_one), " TE entries with single/ambiguous TX_strand)")
