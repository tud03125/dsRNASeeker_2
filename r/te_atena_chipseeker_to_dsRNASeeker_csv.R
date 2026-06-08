#!/usr/bin/env Rscript
# dsRNASeeker TE analysis: standard genomes + custom/T2T-safe annotation
# Corrected 2026-06-06

suppressPackageStartupMessages({
  library(atena)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(regioneR)
  library(ChIPseeker)
  library(BiocParallel)
  library(AnnotationDbi)
  library(S4Vectors)
  library(GenomeInfoDb)
})

args <- commandArgs(trailingOnly=TRUE)

get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1]]
}

as_bool <- function(x, default=FALSE) {
  if (is.null(x)) return(default)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

samplesheet_file <- get_arg("--samplesheet")
case_label       <- get_arg("--case")
control_label    <- get_arg("--control")
outdir           <- get_arg("--outdir")
out_file         <- get_arg("--out")
genome           <- get_arg("--genome")
rmsk_rds         <- get_arg("--rmsk-rds", file.path(outdir, paste0("rmsk_", genome, "_used.rds")))
force_rebuild    <- as_bool(get_arg("--force-rebuild-rmsk", "FALSE"))
single_end       <- !as_bool(get_arg("--paired", "TRUE"))
use_strand       <- as_bool(get_arg("--use-strand", "TRUE"), default=TRUE)
yield_size       <- as.integer(get_arg("--yield-size", "1000000"))
min_max_count    <- as.numeric(get_arg("--min-max-count", "1"))
alpha            <- as.numeric(get_arg("--alpha", "0.10"))
lfc_threshold    <- as.numeric(get_arg("--lfc-threshold", "1"))
shrink_type      <- get_arg("--shrink-type", "ashr")
txdb_package     <- get_arg("--txdb-package", NULL)
orgdb_package    <- get_arg("--orgdb-package", NULL)
txdb_gtf         <- get_arg("--txdb-gtf", NULL)
txdb_rds         <- get_arg("--txdb-rds", file.path(outdir, paste0("TxDb_", genome, ".sqlite")))

if (is.null(samplesheet_file) || is.null(case_label) || is.null(control_label) ||
    is.null(outdir) || is.null(out_file) || is.null(genome)) {
  stop("Required: --samplesheet --case --control --outdir --out --genome")
}

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "counts"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "deseq2"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, "annotation"), recursive=TRUE, showWarnings=FALSE)

message("Advanced TE mode: atena/qtex + DESeq2 + ChIPseeker")
message("Genome: ", genome)
message("Samplesheet: ", samplesheet_file)
message("Output directory: ", outdir)

# -------------------------------------------------------------------------
# TxDb / OrgDb selection
# -------------------------------------------------------------------------
infer_orgdb <- function(genome, orgdb_package=NULL) {
  if (!is.null(orgdb_package) && nzchar(orgdb_package)) return(orgdb_package)
  g <- tolower(genome)
  if (g == "hg38") return("org.Hs.eg.db")
  if (g %in% c("mm10", "mm39", "custom", "t2t", "c57bl_6j_t2t_v1")) return("org.Mm.eg.db")
  stop("Unknown --genome '", genome, "'. Supply --orgdb-package and either --txdb-package or --txdb-gtf.")
}

infer_txdb_package <- function(genome, txdb_package=NULL) {
  if (!is.null(txdb_package) && nzchar(txdb_package)) return(txdb_package)
  g <- tolower(genome)
  if (g == "hg38") return("TxDb.Hsapiens.UCSC.hg38.knownGene")
  if (g == "mm39") return("TxDb.Mmusculus.UCSC.mm39.knownGene")
  if (g == "mm10") return("TxDb.Mmusculus.UCSC.mm10.knownGene")
  stop("No packaged TxDb is defined for --genome '", genome, "'. Supply --txdb-gtf for a custom/T2T assembly.")
}

build_or_load_txdb <- function(txdb_gtf=NULL, txdb_rds=NULL, genome, txdb_package=NULL) {
  use_custom <- !is.null(txdb_gtf) && nzchar(txdb_gtf)
  if (use_custom) {
    if (!file.exists(txdb_gtf)) stop("Custom TxDb GTF does not exist: ", txdb_gtf)
    if (!is.null(txdb_rds) && nzchar(txdb_rds) && file.exists(txdb_rds) && file.info(txdb_rds)$size > 0) {
      message("Loading cached custom TxDb SQLite database: ", txdb_rds)
      txdb_obj <- AnnotationDbi::loadDb(txdb_rds)
      if (!inherits(txdb_obj, "TxDb")) stop("Cached TxDb database did not load as a TxDb object: ", txdb_rds)
      return(txdb_obj)
    }
    message("Building custom TxDb from GTF: ", txdb_gtf)
    organism_name <- if (tolower(genome) %in% c("custom", "t2t", "c57bl_6j_t2t_v1")) "Mus musculus" else NA_character_
    if (requireNamespace("txdbmaker", quietly=TRUE)) {
      txdb_obj <- txdbmaker::makeTxDbFromGFF(file=txdb_gtf, format="gtf", dataSource=txdb_gtf, organism=organism_name)
    } else if (requireNamespace("GenomicFeatures", quietly=TRUE)) {
      txdb_obj <- GenomicFeatures::makeTxDbFromGFF(file=txdb_gtf, format="gtf", dataSource=txdb_gtf, organism=organism_name)
    } else {
      stop("Building a custom TxDb requires either 'txdbmaker' or 'GenomicFeatures'.")
    }
    if (!is.null(txdb_rds) && nzchar(txdb_rds)) {
      dir.create(dirname(txdb_rds), recursive=TRUE, showWarnings=FALSE)
      AnnotationDbi::saveDb(txdb_obj, file=txdb_rds)
      message("Saved custom TxDb SQLite database: ", txdb_rds)
    }
    return(txdb_obj)
  }
  pkg <- infer_txdb_package(genome, txdb_package)
  if (!requireNamespace(pkg, quietly=TRUE)) stop("Required TxDb package not installed: ", pkg, "\nInstall it or supply --txdb-gtf.")
  suppressPackageStartupMessages(library(pkg, character.only=TRUE))
  get(pkg)
}

orgdb_package <- infer_orgdb(genome, orgdb_package)
if (!requireNamespace(orgdb_package, quietly=TRUE)) stop("Required OrgDb package not installed: ", orgdb_package)
suppressPackageStartupMessages(library(orgdb_package, character.only=TRUE))

txdb <- build_or_load_txdb(txdb_gtf=txdb_gtf, txdb_rds=txdb_rds, genome=genome, txdb_package=txdb_package)
annoDb <- orgdb_package
txdb_source_label <- if (!is.null(txdb_gtf) && nzchar(txdb_gtf)) {
  txdb_gtf
} else {
  infer_txdb_package(genome, txdb_package)
}
message("TxDb source: ", txdb_source_label)
message("OrgDb package: ", orgdb_package)

# -------------------------------------------------------------------------
# Samplesheet
# -------------------------------------------------------------------------
samples <- read.delim(samplesheet_file, stringsAsFactors=FALSE, check.names=FALSE)
required_cols <- c("sample_id", "condition", "bam_path")
missing <- setdiff(required_cols, colnames(samples))
if (length(missing) > 0) {
  stop("Samplesheet missing required columns: ", paste(missing, collapse=", "))
}

samples <- samples %>%
  mutate(
    sample_id = as.character(sample_id),
    sample = sample_id,
    condition = as.character(condition),
    group = factor(condition, levels=c(control_label, case_label)),
    bam = bam_path,
    strandedness = if ("strandedness" %in% colnames(.)) strandedness else NA_character_,
    single_end = single_end,
    read_length_bp = NA_integer_
  )

samples <- samples %>% filter(condition %in% c(control_label, case_label))
if (!all(c(control_label, case_label) %in% samples$condition)) {
  stop("Both control and case labels must be present in samplesheet. Found: ",
       paste(unique(samples$condition), collapse=", "))
}

if (any(is.na(samples$group))) {
  stop("Could not set group factor for all samples.")
}

if (!all(file.exists(samples$bam))) {
  missing_bams <- samples$bam[!file.exists(samples$bam)]
  stop("Missing BAM files:\n", paste(missing_bams, collapse="\n"))
}

write.csv(samples, file.path(outdir, paste0("samplesheet_", control_label, "_vs_", case_label, ".csv")), row.names=FALSE)

# -------------------------------------------------------------------------
# RMSK loading/building
# -------------------------------------------------------------------------
validate_rmsk <- function(x) {
  needed_cols <- c("repName", "repClass", "repFamily")
  is_ok <- TRUE
  reasons <- character(0)
  if (!inherits(x, "GRanges")) {
    is_ok <- FALSE
    reasons <- c(reasons, "Object is not a GRanges.")
  } else {
    if (length(x) == 0) {
      is_ok <- FALSE
      reasons <- c(reasons, "GRanges is empty.")
    }
    missing_cols <- setdiff(needed_cols, names(S4Vectors::mcols(x)))
    if (length(missing_cols) > 0) {
      is_ok <- FALSE
      reasons <- c(reasons, paste0("Missing metadata columns: ", paste(missing_cols, collapse=", ")))
    }
  }
  list(ok=is_ok, reasons=reasons)
}

load_or_build_rmsk <- function(rds_path, genome, force_rebuild=FALSE) {
  if (!force_rebuild && file.exists(rds_path)) {
    message("Loading cached RMSK: ", rds_path)
    x <- readRDS(rds_path)
    chk <- validate_rmsk(x)
    if (chk$ok) return(x)
    message("Cached RMSK failed validation; rebuilding.")
    for (rr in chk$reasons) message("  - ", rr)
  }

  if (tolower(genome) %in% c("custom", "t2t", "c57bl_6j_t2t_v1")) {
    stop("For a custom/T2T assembly, --rmsk-rds must point to an existing validated GRanges RDS with repName, repClass, and repFamily metadata.")
  }
  message("Building RMSK annotation with atena::annotaTEs(genome='", genome, "') ...")
  x <- annotaTEs(genome=genome, parsefun=rmskbasicparser)
  saveRDS(x, rds_path)
  message("Saved RMSK RDS: ", rds_path)
  x
}

rmsk <- load_or_build_rmsk(rmsk_rds, genome, force_rebuild=force_rebuild)

rmsk_qc <- tibble(
  metric = c("class", "length", "seqlevels_n", "metadata_cols"),
  value = c(
    class(rmsk)[1],
    as.character(length(rmsk)),
    as.character(length(GenomeInfoDb::seqlevels(rmsk))),
    paste(names(S4Vectors::mcols(rmsk)), collapse=";")
  )
)
write.csv(rmsk_qc, file.path(outdir, paste0("rmsk_", genome, "_qc.csv")), row.names=FALSE)

rmsk_seq <- GenomeInfoDb::seqlevels(rmsk)
txdb_seq <- GenomeInfoDb::seqlevels(txdb)
shared_seq <- intersect(rmsk_seq, txdb_seq)
write.csv(data.frame(metric=c("rmsk_seqlevels_n","txdb_seqlevels_n","shared_seqlevels_n","shared_seqlevels"), value=c(length(rmsk_seq),length(txdb_seq),length(shared_seq),paste(shared_seq,collapse=";"))), file.path(outdir,paste0("coordinate_compatibility_",genome,".csv")), row.names=FALSE)
if (length(shared_seq) == 0) stop("RepeatMasker GRanges and TxDb have zero shared sequence names; FASTA/GTF/RMSK are not coordinate-compatible.")
message("Shared RMSK/TxDb seqlevels: ", length(shared_seq))

# -------------------------------------------------------------------------
# atena / qtex quantification
# -------------------------------------------------------------------------
message("Creating TEtranscriptsParam ...")
ttpar <- TEtranscriptsParam(
  bfl = samples$bam,
  teFeatures = rmsk,
  singleEnd = single_end,
  ignoreStrand = !use_strand
)

pheno <- as.data.frame(samples[, c("sample", "condition", "group", "strandedness", "single_end", "read_length_bp")])
rownames(pheno) <- samples$sample

qtex_rds <- file.path(outdir, "TE_qtex_summarized_experiment.rds")
counts_csv <- file.path(outdir, "counts", "TE_raw_counts_all_samples.csv")

if (file.exists(qtex_rds) && file.info(qtex_rds)$size > 0) {
  message("Loading cached qtex SummarizedExperiment: ", qtex_rds)
  emq <- readRDS(qtex_rds)

  # Re-write counts CSV if missing, without rerunning qtex().
  if (!file.exists(counts_csv) || file.info(counts_csv)$size == 0) {
    message("Writing TE count matrix from cached qtex object: ", counts_csv)
    write.csv(as.data.frame(assay(emq)), counts_csv)
  }
} else {
  # -------------------------------------------------------------------------
  # atena custom-assembly seqlevel patch
  #
  # atena normally calls GenomeInfoDb::seqlevelsStyle() to reconcile BAM,
  # TE, and gene chromosome naming conventions. That inference fails for
  # valid custom-accession contigs such as CAXLPS010000001.1 even when the
  # BAM and annotations already use exactly matching sequence names.
  #
  # For custom/T2T assemblies, compare exact seqlevel names first. Style
  # conversion is attempted only when exact matching is insufficient.
  # -------------------------------------------------------------------------
  
  if (tolower(genome) %in% c(
      "custom",
      "t2t",
      "c57bl_6j_t2t_v1"
  )) {
      message(
          "Installing targeted atena custom-assembly seqlevel patch"
      )
  
      patched_joinTEsGenes <- function(
          teFeatures,
          geneFeatures,
          geneFeaturesobjname = deparse(substitute(geneFeatures))
      ) {
          if (
              !methods::is(geneFeatures, "GRanges") &&
              !methods::is(geneFeatures, "GRangesList")
          ) {
              stop(
                  sprintf(
                      paste0(
                          "gene features object '%s' should be either ",
                          "a 'GRanges' or a 'GRangesList' object."
                      ),
                      geneFeaturesobjname
                  )
              )
          }
  
          if (is.null(names(geneFeatures))) {
              stop(
                  sprintf(
                      "gene features object '%s' has no names().",
                      geneFeaturesobjname
                  )
              )
          }
  
          if (any(names(geneFeatures) %in% names(teFeatures))) {
              stop(
                  "Gene features have identifiers in common with ",
                  "the TE features."
              )
          }
  
          if (length(geneFeatures) == 0L) {
              stop(
                  sprintf(
                      "gene features object '%s' is empty.",
                      geneFeaturesobjname
                  )
              )
          }
  
          te_levels <- GenomeInfoDb::seqlevels(teFeatures)
          gene_levels <- GenomeInfoDb::seqlevels(geneFeatures)
          exact_shared <- intersect(te_levels, gene_levels)
  
          if (length(exact_shared) > 0L) {
              message(
                  "atena .joinTEsGenes: exact TE/gene seqlevel names ",
                  "already overlap; skipping seqlevelsStyle conversion"
              )
          } else {
              te_style <- tryCatch(
                  GenomeInfoDb::seqlevelsStyle(teFeatures),
                  error = function(e) character(0)
              )
  
              if (length(te_style) == 0L) {
                  stop(
                      "TE and gene features have no exact shared sequence ",
                      "names, and the custom sequence naming style cannot ",
                      "be inferred."
                  )
              }
  
              GenomeInfoDb::seqlevelsStyle(geneFeatures) <-
                  te_style[[1]]
          }
  
          combined_levels <- unique(
              c(
                  GenomeInfoDb::seqlevels(teFeatures),
                  GenomeInfoDb::seqlevels(geneFeatures)
              )
          )
  
          GenomeInfoDb::seqlevels(teFeatures) <- combined_levels
          GenomeInfoDb::seqlevels(geneFeatures) <- combined_levels
  
          features <- c(teFeatures, geneFeatures)
  
          te_mask <- S4Vectors::Rle(
              rep(
                  FALSE,
                  length(teFeatures) + length(geneFeatures)
              )
          )
  
          te_mask[seq_along(teFeatures)] <- TRUE
  
          S4Vectors::mcols(features)$isTE <- te_mask
  
          features
      }
  
      patched_matchSeqinfo <- function(
          gal,
          features,
          verbose = TRUE
      ) {
          # Use formal S4 inheritance checks rather than exact class-name
          # membership. qtex supplies subclasses such as
          # CompressedGRangesList, which inherits from GRangesList.
          valid_alignment <- (
              methods::is(gal, "GAlignments") ||
              methods::is(gal, "GAlignmentPairs") ||
              methods::is(gal, "GAlignmentsList")
          )
  
          valid_features <- (
              methods::is(features, "GRanges") ||
              methods::is(features, "GRangesList")
          )
  
          if (!valid_alignment || !valid_features) {
              stop(
                  "Unexpected object classes in atena .matchSeqinfo: ",
                  "gal=", paste(class(gal), collapse = "/"),
                  "; features=", paste(class(features), collapse = "/"),
                  ". S4 inheritance checks failed."
              )
          }
  
          gal_levels <- GenomeInfoDb::seqlevels(gal)
          feature_levels <- GenomeInfoDb::seqlevels(features)
  
          common_chr <- intersect(
              gal_levels,
              feature_levels
          )
  
          if (length(common_chr) > 0L) {
              if (verbose && !isTRUE(getOption("dsRNASeeker.atena_match_message_shown"))) {
                  message(
                      "atena .matchSeqinfo: found ",
                      length(common_chr),
                      " exact BAM/annotation seqlevel matches; ",
                      "skipping seqlevelsStyle conversion"
                  )
                  options(dsRNASeeker.atena_match_message_shown=TRUE)
              }
          } else {
              gal_style <- tryCatch(
                  GenomeInfoDb::seqlevelsStyle(gal),
                  error = function(e) character(0)
              )
  
              feature_style <- tryCatch(
                  GenomeInfoDb::seqlevelsStyle(features),
                  error = function(e) character(0)
              )
  
              if (
                  length(gal_style) > 0L &&
                  length(feature_style) > 0L &&
                  length(intersect(gal_style, feature_style)) > 0L
              ) {
                  return(gal)
              }
  
              if (length(feature_style) == 0L) {
                  stop(
                      "The BAM and feature annotations have no exact shared ",
                      "sequence names, and the custom annotation style cannot ",
                      "be inferred."
                  )
              }
  
              GenomeInfoDb::seqlevelsStyle(gal) <-
                  feature_style[[1]]
  
              gal_levels <- GenomeInfoDb::seqlevels(gal)
              feature_levels <- GenomeInfoDb::seqlevels(features)
  
              common_chr <- intersect(
                  gal_levels,
                  feature_levels
              )
  
              if (length(common_chr) == 0L) {
                  stop(
                      "No shared sequence names remain between the BAM ",
                      "and annotation after attempted style conversion."
                  )
              }
          }
  
          gal_lengths <- GenomeInfoDb::seqlengths(gal)[common_chr]
          feature_lengths <-
              GenomeInfoDb::seqlengths(features)[common_chr]
  
          comparable <- (
              !is.na(gal_lengths) &
              !is.na(feature_lengths)
          )
  
          length_mismatch <- (
              comparable &
              gal_lengths != feature_lengths
          )
  
          if (any(length_mismatch)) {
              mismatched_chr <- common_chr[length_mismatch]
  
              if (verbose) {
                  message(
                      "Sequence-length mismatch between BAM and annotation ",
                      "for: ",
                      paste(mismatched_chr, collapse = ", "),
                      ". These sequence levels will be discarded."
                  )
              }
  
              common_chr <- common_chr[!length_mismatch]
          }
  
          if (length(common_chr) == 0L) {
              stop(
                  "None of the shared BAM/annotation sequence levels ",
                  "have compatible sequence lengths."
              )
          }
  
          gal <- GenomeInfoDb::keepSeqlevels(
              gal,
              common_chr,
              pruning.mode = "coarse"
          )
  
          feature_seqinfo <- GenomeInfoDb::seqinfo(features)[common_chr]
  
          old_order <- match(
              common_chr,
              GenomeInfoDb::seqlevels(gal)
          )
  
          GenomeInfoDb::seqinfo(
              gal,
              new2old = old_order,
              pruning.mode = "coarse"
          ) <- feature_seqinfo
  
          gal
      }
  
      assignInNamespace(
          ".joinTEsGenes",
          patched_joinTEsGenes,
          ns = "atena"
      )
  
      assignInNamespace(
          ".matchSeqinfo",
          patched_matchSeqinfo,
          ns = "atena"
      )
  
      message(
          "Targeted atena custom-assembly seqlevel patch installed"
      )
  }
  message("Running atena::qtex() ...")
  emq <- qtex(
    ttpar,
    phenodata = pheno,
    mode = ovUnion,
    yieldSize = yield_size,
    BPPARAM = SerialParam(progressbar=FALSE)
  )

  # Try to force expected sample names when qtex preserves BAM basenames instead.
  if (ncol(assay(emq)) == nrow(samples)) {
    colnames(assay(emq)) <- samples$sample
  }

  saveRDS(emq, qtex_rds)
  message("Saved qtex SummarizedExperiment: ", qtex_rds)

  write.csv(as.data.frame(assay(emq)), counts_csv)
}

# -------------------------------------------------------------------------
# Differential analysis
# -------------------------------------------------------------------------
run_deseq_te <- function(emq_obj, samples_df, condition1, condition2,
                         alpha=0.10, lfc_threshold=1, shrink_type="ashr",
                         min_max_count=1) {
  countData <- assay(emq_obj)

  sample_names <- samples_df$sample
  if (!all(sample_names %in% colnames(countData))) {
    stop("Not all sample names are present in qtex count matrix. Missing: ",
         paste(setdiff(sample_names, colnames(countData)), collapse=", "),
         "\nCount matrix columns: ", paste(colnames(countData), collapse=", "))
  }

  colDt <- as.data.frame(samples_df)
  rownames(colDt) <- colDt$sample
  colDt <- colDt[sample_names, , drop=FALSE]

  subset_countData <- countData[, sample_names, drop=FALSE]
  subset_countData <- subset_countData[apply(subset_countData, 1, max, na.rm=TRUE) > min_max_count, , drop=FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = round(subset_countData),
    colData = colDt,
    design = ~ group
  )
  dds$group <- relevel(dds$group, ref=condition1)
  dds <- DESeq(dds, quiet=TRUE)

  res <- results(dds, contrast=c("group", condition2, condition1), alpha=alpha)
  norm_counts <- counts(dds, normalized=TRUE)

  mean_group1 <- rowMeans(norm_counts[, colData(dds)$group == condition1, drop=FALSE])
  mean_group2 <- rowMeans(norm_counts[, colData(dds)$group == condition2, drop=FALSE])

  res_df <- as.data.frame(res) %>%
    rownames_to_column("TE_id") %>%
    mutate(mean_group1 = mean_group1[TE_id],
           mean_group2 = mean_group2[TE_id])

  shrunk <- tryCatch({
    if (requireNamespace("ashr", quietly=TRUE) && shrink_type == "ashr") {
      lfcShrink(dds,
                contrast=c("group", condition2, condition1),
                res=res,
                alpha=alpha,
                lfcThreshold=lfc_threshold,
                type="ashr")
    } else {
      message("ashr not available or not requested; using unshrunk DESeq2 results.")
      res
    }
  }, error=function(e) {
    message("lfcShrink failed; using unshrunk results. Reason: ", conditionMessage(e))
    res
  })

  shrunk_df <- as.data.frame(shrunk) %>%
    rownames_to_column("TE_id") %>%
    mutate(mean_group1 = mean_group1[TE_id],
           mean_group2 = mean_group2[TE_id])

  list(
    dds=dds,
    res=res_df,
    shrunk=shrunk_df,
    norm_counts=as.data.frame(norm_counts) %>% rownames_to_column("TE_id")
  )
}

message("Running DESeq2: ", control_label, " vs ", case_label)
de_out <- run_deseq_te(
  emq_obj=emq,
  samples_df=samples,
  condition1=control_label,
  condition2=case_label,
  alpha=alpha,
  lfc_threshold=lfc_threshold,
  shrink_type=shrink_type,
  min_max_count=min_max_count
)

comparison <- paste0(control_label, "_vs_", case_label)

write.csv(de_out$res, file.path(outdir, "deseq2", paste0("DESeq2_results_unshrunk_", comparison, ".csv")), row.names=FALSE)
write.csv(de_out$shrunk, file.path(outdir, "deseq2", paste0("DESeq2_results_shrunk_", comparison, ".csv")), row.names=FALSE)
write.csv(de_out$norm_counts, file.path(outdir, "deseq2", paste0("DESeq2_normalized_counts_", comparison, ".csv")), row.names=FALSE)

sig_all <- de_out$shrunk %>%
  filter(!is.na(padj), padj < alpha, abs(log2FoldChange) > lfc_threshold)

sig_up <- de_out$shrunk %>%
  filter(!is.na(padj), padj < alpha, log2FoldChange > lfc_threshold)

# -------------------------------------------------------------------------
# RepeatMasker metadata extraction
# -------------------------------------------------------------------------
extract_rmsk_annotation <- function(rmsk_obj, te_ids) {
  nm <- names(ranges(rmsk_obj))
  if (is.null(nm)) nm <- names(rmsk_obj)
  rr <- rmsk_obj[nm %in% te_ids]
  ann <- as.data.frame(rr)

  if (nrow(ann) == 0) {
    return(data.frame(
      TE_id=character(), seqnames=character(), start=integer(), end=integer(),
      width=integer(), strand=character(), repName=character(),
      repClass=character(), repFamily=character(), repStart=integer(),
      repEnd=integer(), stringsAsFactors=FALSE
    ))
  }

  if (!"repStart" %in% colnames(ann)) ann$repStart <- NA
  if (!"repEnd" %in% colnames(ann)) ann$repEnd <- NA
  if (!"strand" %in% colnames(ann)) ann$strand <- as.character(strand(rr))
  ann$TE_id <- rownames(ann)

  keep <- intersect(
    c("TE_id", "seqnames", "start", "end", "width", "strand",
      "repName", "repClass", "repFamily", "repStart", "repEnd"),
    colnames(ann)
  )
  ann[, keep, drop=FALSE]
}

rmsk_ann_all <- extract_rmsk_annotation(rmsk, sig_all$TE_id)
rmsk_ann_up  <- extract_rmsk_annotation(rmsk, sig_up$TE_id)

sig_all_annot <- sig_all %>%
  left_join(rmsk_ann_all, by="TE_id") %>%
  arrange(desc(log2FoldChange))

sig_up_annot <- sig_up %>%
  left_join(rmsk_ann_up, by="TE_id") %>%
  arrange(desc(log2FoldChange))

write.csv(sig_all_annot, file.path(outdir, "annotation", paste0("TEtranscript_nn_", comparison, "_results.csv")), row.names=FALSE)
write.csv(sig_up_annot, file.path(outdir, "annotation", paste0("TEtranscript_", comparison, "_results.csv")), row.names=FALSE)

# -------------------------------------------------------------------------
# ChIPseeker-style genomic annotation
# -------------------------------------------------------------------------
is_custom_genome <- tolower(genome) %in% c(
  "custom", "t2t", "c57bl_6j_t2t_v1"
)

empty_peak_annotation <- function(te_ids, annotation_text=NA_character_) {
  te_ids <- as.character(te_ids)
  n <- length(te_ids)

  data.frame(
    TE_id = te_ids,
    annotation = rep(annotation_text, n),
    geneChr = rep(NA_character_, n),
    geneStart = rep(NA_integer_, n),
    geneEnd = rep(NA_integer_, n),
    geneLength = rep(NA_integer_, n),
    geneStrand = rep(NA_character_, n),
    geneId = rep(NA_character_, n),
    transcriptId = rep(NA_character_, n),
    distanceToTSS = rep(NA_integer_, n),
    SYMBOL = rep(NA_character_, n),
    GENENAME = rep(NA_character_, n),
    ENSEMBL = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

make_peak_granges <- function(te_table) {
  required <- c("seqnames", "start", "end", "TE_id", "strand")
  missing_cols <- setdiff(required, colnames(te_table))
  if (length(missing_cols) > 0L) {
    stop(
      "TE table is missing columns required for genomic annotation: ",
      paste(missing_cols, collapse=", ")
    )
  }

  x <- te_table[, required, drop=FALSE]
  x$seqnames <- as.character(x$seqnames)
  x$start <- suppressWarnings(as.integer(x$start))
  x$end <- suppressWarnings(as.integer(x$end))
  x$TE_id <- as.character(x$TE_id)
  x$strand <- as.character(x$strand)
  x$strand[!x$strand %in% c("+", "-", "*")] <- "*"

  valid <- (
    !is.na(x$seqnames) & nzchar(x$seqnames) &
    !is.na(x$start) & !is.na(x$end) &
    x$start >= 1L & x$end >= x$start &
    !is.na(x$TE_id) & nzchar(x$TE_id)
  )
  x <- x[valid, , drop=FALSE]

  if (nrow(x) == 0L) return(GenomicRanges::GRanges())

  gr <- GenomicRanges::GRanges(
    seqnames = x$seqnames,
    ranges = IRanges::IRanges(start=x$start, end=x$end),
    strand = x$strand
  )
  names(gr) <- make.unique(x$TE_id)
  S4Vectors::mcols(gr)$TE_id <- x$TE_id
  gr
}

annotate_te_regions <- function(te_table, out_prefix, txdb, annoDb,
                                custom_genome=FALSE) {
  output_csv <- file.path(
    outdir, "annotation", paste0(out_prefix, "_peak_annotation.csv")
  )

  if (nrow(te_table) == 0L) {
    message("No TE rows to annotate for ", out_prefix)
    return(NULL)
  }

  # Standard reference genomes deliberately follow the original, validated
  # hg38/mm10/mm39 ChIPseeker route. No T2T filtering or placeholder table
  # is applied to these assemblies.
  if (!custom_genome) {
    required <- c("seqnames", "start", "end", "TE_id", "strand", "repName")
    missing_cols <- setdiff(required, colnames(te_table))
    if (length(missing_cols) > 0L) {
      stop(
        "TE table is missing columns required for ChIPseeker annotation: ",
        paste(missing_cols, collapse=", ")
      )
    }

    gr_input <- te_table %>%
      dplyr::select(
        seqnames, start, end, TE_id, strand, repName
      ) %>%
      dplyr::filter(
        !is.na(seqnames),
        !is.na(start),
        !is.na(end),
        !is.na(TE_id)
      )

    if (nrow(gr_input) == 0L) {
      message("No valid genomic ranges to annotate for ", out_prefix)
      return(NULL)
    }

    gr <- regioneR::toGRanges(gr_input)

    peak_anno <- ChIPseeker::annotatePeak(
      peak = gr,
      TxDb = txdb,
      annoDb = annoDb,
      verbose = FALSE
    )

    anno_df <- as.data.frame(peak_anno)
    anno_df <- anno_df[, !duplicated(colnames(anno_df)), drop=FALSE]

    if (!"TE_id" %in% colnames(anno_df)) {
      candidate_cols <- intersect(
        c("name", "peak", "V4"),
        colnames(anno_df)
      )
      if (length(candidate_cols) > 0L) {
        anno_df$TE_id <- as.character(
          anno_df[[candidate_cols[[1]]]]
        )
      } else if (nrow(anno_df) == nrow(gr_input)) {
        anno_df$TE_id <- as.character(gr_input$TE_id)
      } else {
        stop(
          "Standard-genome ChIPseeker output lacks a recoverable TE_id."
        )
      }
    }

    write.csv(anno_df, output_csv, row.names=FALSE)
    return(anno_df)
  }

  # Custom/T2T route. These assemblies can contain repeat-bearing contigs
  # absent from the supplied gene GTF/TxDb. Only TxDb-supported contigs are
  # sent to ChIPseeker; all other significant TEs remain in the output with
  # an explicit annotation status.
  gr_all <- make_peak_granges(te_table)
  if (length(gr_all) == 0L) {
    message("No valid genomic ranges to annotate for ", out_prefix)
    return(NULL)
  }

  txdb_levels <- GenomeInfoDb::seqlevels(txdb)
  peak_levels <- as.character(GenomeInfoDb::seqnames(gr_all))
  supported <- peak_levels %in% txdb_levels

  unsupported_ids <- S4Vectors::mcols(gr_all)$TE_id[!supported]
  unsupported_df <- empty_peak_annotation(
    unsupported_ids,
    "Unannotated: contig absent from TxDb"
  )

  gr <- gr_all[supported]

  if (length(gr) == 0L) {
    message(
      "No TE ranges for ", out_prefix,
      " occur on contigs represented in the TxDb; writing unannotated rows."
    )
    write.csv(unsupported_df, output_csv, row.names=FALSE)
    return(unsupported_df)
  }

  shared_levels <- intersect(
    unique(as.character(GenomeInfoDb::seqnames(gr))),
    txdb_levels
  )

  gr <- GenomeInfoDb::keepSeqlevels(
    gr,
    shared_levels,
    pruning.mode="coarse"
  )

  txdb_si <- GenomeInfoDb::seqinfo(txdb)[shared_levels]
  GenomeInfoDb::seqinfo(gr) <- txdb_si

  peak_anno <- tryCatch(
    ChIPseeker::annotatePeak(
      peak=gr,
      TxDb=txdb,
      annoDb=NULL,
      verbose=FALSE
    ),
    error=function(e) {
      message(
        "ChIPseeker annotation failed for ", out_prefix,
        "; retaining TE rows with NA genomic annotation. Reason: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(peak_anno)) {
    supported_df <- empty_peak_annotation(
      S4Vectors::mcols(gr)$TE_id,
      "Unannotated: ChIPseeker failed"
    )
  } else {
    supported_df <- as.data.frame(peak_anno)
    supported_df <- supported_df[
      , !duplicated(colnames(supported_df)),
      drop=FALSE
    ]

    if (!"TE_id" %in% colnames(supported_df)) {
      candidate_cols <- intersect(
        c("name", "peak", "V4"),
        colnames(supported_df)
      )

      if (length(candidate_cols) > 0L) {
        supported_df$TE_id <- as.character(
          supported_df[[candidate_cols[[1]]]]
        )
      } else if (nrow(supported_df) == length(gr)) {
        supported_df$TE_id <- as.character(
          S4Vectors::mcols(gr)$TE_id
        )
      } else {
        stop(
          "Custom-genome ChIPseeker output lacks a recoverable TE_id."
        )
      }
    }
  }

  anno_df <- dplyr::bind_rows(
    supported_df,
    unsupported_df
  )

  if (nrow(anno_df) > 0L) {
    anno_df <- anno_df[
      !duplicated(anno_df$TE_id),
      ,
      drop=FALSE
    ]
  }

  write.csv(anno_df, output_csv, row.names=FALSE)
  anno_df
}

message("Annotating genomic context with ChIPseeker ...")
all_peak_anno <- annotate_te_regions(
  sig_all_annot,
  paste0(comparison, "_all_sig"),
  txdb,
  annoDb,
  custom_genome=is_custom_genome
)
up_peak_anno <- annotate_te_regions(
  sig_up_annot,
  paste0(comparison, "_up_sig"),
  txdb,
  annoDb,
  custom_genome=is_custom_genome
)

safe_merge_peak <- function(expr_df, peak_df, outfile) {
  if (is.null(peak_df) || nrow(expr_df) == 0) {
    write.csv(expr_df, outfile, row.names=FALSE)
    return(expr_df)
  }

  if (!"TE_id" %in% colnames(peak_df)) {
    message("Peak annotation lacks TE_id; writing expression table without peak merge: ", outfile)
    write.csv(expr_df, outfile, row.names=FALSE)
    return(expr_df)
  }

  peak_df$TE_id <- as.character(peak_df$TE_id)
  expr_df$TE_id <- as.character(expr_df$TE_id)

  keep_cols <- intersect(
    c("TE_id", "annotation", "geneChr", "geneStart", "geneEnd", "geneLength",
      "geneStrand", "geneId", "transcriptId", "distanceToTSS", "SYMBOL",
      "GENENAME", "ENSEMBL"),
    colnames(peak_df)
  )

  merged <- expr_df %>%
    left_join(peak_df[, keep_cols, drop=FALSE], by="TE_id")

  write.csv(merged, outfile, row.names=FALSE)
  merged
}

merged_all <- safe_merge_peak(
  sig_all_annot,
  all_peak_anno,
  file.path(outdir, "annotation", paste0("TE_expression_annotation_", comparison, "_all_sig.csv"))
)

merged_up <- safe_merge_peak(
  sig_up_annot,
  up_peak_anno,
  file.path(outdir, "annotation", paste0("TE_expression_annotation_", comparison, "_up_sig.csv"))
)

# dsRNASeeker-facing canonical CSV. This is the file consumed by the Python pipeline.
write.csv(merged_all, out_file, row.names=FALSE)

# Intron subset
extract_introns <- function(expr_df, peak_df, outfile) {
  if (is.null(peak_df) || nrow(expr_df) == 0) {
    write.csv(expr_df[0, , drop=FALSE], outfile, row.names=FALSE)
    return(invisible(NULL))
  }
  if (!"TE_id" %in% colnames(peak_df) || !"annotation" %in% colnames(peak_df)) {
    write.csv(expr_df[0, , drop=FALSE], outfile, row.names=FALSE)
    return(invisible(NULL))
  }
  peak_df$TE_id <- as.character(peak_df$TE_id)
  expr_df$TE_id <- as.character(expr_df$TE_id)

  intron_hits <- peak_df %>%
    dplyr::filter(grepl("Intron", annotation, ignore.case=TRUE)) %>%
    dplyr::select(TE_id) %>%
    dplyr::distinct()

  intron_df <- expr_df %>% inner_join(intron_hits, by="TE_id")
  write.csv(intron_df, outfile, row.names=FALSE)
  invisible(intron_df)
}

extract_introns(merged_all, all_peak_anno, file.path(outdir, "annotation", paste0("intron_TE_", comparison, "_all_sig.csv")))
extract_introns(merged_up,  up_peak_anno,  file.path(outdir, "annotation", paste0("intron_TE_", comparison, "_up_sig.csv")))

summary_df <- tibble(
  metric = c("n_total_TE_features", "n_tested_after_count_filter", "n_sig_all_directions", "n_sig_up_only"),
  value = c(nrow(assay(emq)), nrow(de_out$shrunk), nrow(sig_all_annot), nrow(sig_up_annot))
)
write.csv(summary_df, file.path(outdir, paste0("run_summary_", comparison, ".csv")), row.names=FALSE)

message("Done advanced TE analysis.")
message("Canonical dsRNASeeker TE CSV: ", out_file)
