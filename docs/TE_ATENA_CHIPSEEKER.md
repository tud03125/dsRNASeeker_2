# Advanced TE module: atena + DESeq2 + ChIPseeker

This patch restores the TE-analysis strategy used in the user's original scripts:
  BAMs -> atena::qtex/TEtranscriptsParam -> DESeq2 -> ChIPseeker annotation -> dsRNASeeker CSV

## Files to replace/add

Replace:
  modules/te_analysis.py
  main.py
  environment_dsRNASeeker_2.yml

Add:
  r/te_atena_chipseeker_to_dsRNASeeker_csv.R

Keep:
  r/te_deseq2_to_dsRNASeeker_csv.R
as fallback for --te-mode simple.

## New workflow options

  --te-mode advanced              default; atena + DESeq2 + ChIPseeker
  --te-mode simple                old fallback; featureCounts + DESeq2
  --te-genome hg38|mm39|mm10      required for advanced mode
  --te-rmsk-rds PATH              optional cached RMSK GRanges RDS
  --te-force-rebuild-rmsk         rebuild RMSK with atena::annotaTEs()
  --te-use-strand                 default TRUE
  --te-ignore-strand              use ignoreStrand=TRUE inside atena
  --te-yield-size 1000000
  --te-min-max-count 1
  --te-shrink-type ashr|none
  --te-txdb-package PACKAGE       optional override
  --te-orgdb-package PACKAGE      optional override

## GSE59717 example

Add these to the workflow CLI:

  --te-mode advanced \
  --te-genome hg38 \
  --te-rmsk-rds /rs01/home/levinm/rmsk_hg38_used.rds \
  --te-use-strand \
  --te-padj-max 0.05 \
  --te-lfc-min 1

## GSE308489 mm39 example

  --te-mode advanced \
  --te-genome mm39 \
  --te-rmsk-rds /rs01/home/levinm/rmsk_mm39_used.rds \
  --te-use-strand \
  --te-padj-max 0.05 \
  --te-lfc-min 1

## Output contract

The advanced R script writes the canonical dsRNASeeker input file:

  02_te/annotation/TE_expression_annotation_<CONTROL>_vs_<CASE>_all_sig.dsRNASeeker.csv

It also writes companion outputs analogous to the original per-dataset scripts:

  TE_qtex_summarized_experiment.rds
  counts/TE_raw_counts_all_samples.csv
  deseq2/DESeq2_results_unshrunk_<CONTROL>_vs_<CASE>.csv
  deseq2/DESeq2_results_shrunk_<CONTROL>_vs_<CASE>.csv
  annotation/<CONTROL>_vs_<CASE>_all_sig_peak_annotation.csv
  annotation/TE_expression_annotation_<CONTROL>_vs_<CASE>_all_sig.csv
  annotation/intron_TE_<CONTROL>_vs_<CASE>_all_sig.csv

## Installation

The environment YAML now includes the Bioconductor packages required for advanced mode.
If conda cannot solve all Bioconductor packages on FCCC, install missing R packages
inside the active dsRNASeeker_2 environment with BiocManager, for example:

  Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")'
  Rscript -e 'BiocManager::install(c("atena","ChIPseeker","regioneR","TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm39.knownGene","org.Hs.eg.db","org.Mm.eg.db","ashr"), ask=FALSE, update=FALSE)'
