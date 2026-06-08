from __future__ import annotations

from pathlib import Path
import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _featurecounts_strand_code(strandedness: str) -> str:
    s = (strandedness or "auto").lower()
    if s in {"reverse", "fr-firststrand", "firststrand"}:
        return "2"
    if s in {"forward", "fr-secondstrand", "secondstrand"}:
        return "1"
    return "0"


def _bool_flag(x: bool) -> str:
    return "TRUE" if bool(x) else "FALSE"


def run_te_analysis(args, bam_samplesheet: str | Path, te_outdir: str | Path) -> Path:
    """Run the workflow TE module.

    Two modes are supported:
      --te-mode advanced  : atena/qtex + DESeq2 + ChIPseeker annotation.
                            This is the mode that mirrors the user's historical
                            GSE59717/GSE308489 TE-analysis scripts.
      --te-mode simple    : featureCounts + DESeq2 + simple TE-GTF metadata.
                            This remains as a lightweight fallback.

    Both modes produce the same final dsRNASeeker-facing CSV:
      02_te/annotation/TE_expression_annotation_<CONTROL>_vs_<CASE>_all_sig.dsRNASeeker.csv
    """
    te_outdir = ensure_dir(te_outdir)
    annot_dir = ensure_dir(te_outdir / "annotation")
    final_csv = annot_dir / f"TE_expression_annotation_{args.control_label}_vs_{args.case_label}_all_sig.dsRNASeeker.csv"

    if is_nonempty_file(final_csv) and not getattr(args, "force", False):
        step(f"Step 2/6 TE analysis: reusing {final_csv}")
        return final_csv

    mode = getattr(args, "te_mode", "advanced")
    if mode == "advanced":
        return run_te_analysis_atena_chipseeker(args, bam_samplesheet, te_outdir, final_csv)
    if mode == "simple":
        return run_te_analysis_featurecounts(args, bam_samplesheet, te_outdir, final_csv)

    raise ValueError(f"Unknown --te-mode: {mode}")


def run_te_analysis_atena_chipseeker(args, bam_samplesheet: str | Path, te_outdir: Path, final_csv: Path) -> Path:
    """Run atena/qtex + DESeq2 + ChIPseeker TE annotation."""
    ensure_dir(te_outdir / "counts")
    ensure_dir(te_outdir / "deseq2")
    ensure_dir(te_outdir / "annotation")

    sample_meta = te_outdir / "samplesheet.bam.tsv"
    df = pd.read_csv(bam_samplesheet, sep="\t")
    df.to_csv(sample_meta, sep="\t", index=False)

    rscript = Path(__file__).resolve().parents[1] / "r" / "te_atena_chipseeker_to_dsRNASeeker_csv.R"
    log = Path(args.output_dir) / "pipeline_info" / "logs" / "TE_atena_ChIPseeker.log"

    genome = getattr(args, "te_genome", None)
    if not genome:
        raise ValueError(
            "--te-mode advanced requires --te-genome "
            "(for example hg38, mm39, mm10, or custom)."
        )

    genome_key = str(genome).strip().lower()
    custom_genomes = {"custom", "t2t", "c57bl_6j_t2t_v1"}

    # A custom/T2T assembly must never fall back to atena's built-in UCSC
    # annotations. Both the repeat annotation and transcript annotation must
    # come from the same coordinate system as the supplied FASTA/BAMs.
    if genome_key in custom_genomes:
        if not getattr(args, "te_rmsk_rds", None):
            raise ValueError(
                "Custom/T2T advanced TE mode requires --te-rmsk-rds pointing "
                "to an assembly-matched GRanges RDS."
            )
        if not getattr(args, "te_txdb_gtf", None):
            raise ValueError(
                "Custom/T2T advanced TE mode requires --te-txdb-gtf pointing "
                "to an assembly-matched gene GTF."
            )
        if getattr(args, "te_force_rebuild_rmsk", False):
            raise ValueError(
                "--te-force-rebuild-rmsk is not supported for custom/T2T "
                "assemblies. Build the custom RMSK RDS explicitly and pass it "
                "with --te-rmsk-rds."
            )

    cmd = [
        args.rscript_exe, str(rscript),
        "--samplesheet", str(sample_meta),
        "--case", args.case_label,
        "--control", args.control_label,
        "--outdir", str(te_outdir),
        "--out", str(final_csv),
        "--genome", str(genome),
        "--paired", _bool_flag(getattr(args, "paired", True)),
        "--use-strand", _bool_flag(getattr(args, "te_use_strand", True)),
        "--yield-size", str(getattr(args, "te_yield_size", 1000000)),
        "--min-max-count", str(getattr(args, "te_min_max_count", 1)),
        "--alpha", str(getattr(args, "te_padj_max", 0.10)),
        "--lfc-threshold", str(getattr(args, "te_lfc_min", 1.0)),
        "--shrink-type", str(getattr(args, "te_shrink_type", "ashr")),
    ]

    if getattr(args, "te_rmsk_rds", None):
        cmd += ["--rmsk-rds", str(args.te_rmsk_rds)]
    if getattr(args, "te_force_rebuild_rmsk", False):
        cmd += ["--force-rebuild-rmsk", "TRUE"]
    if getattr(args, "te_txdb_package", None):
        cmd += ["--txdb-package", str(args.te_txdb_package)]
    if getattr(args, "te_txdb_gtf", None):
        cmd += ["--txdb-gtf", str(args.te_txdb_gtf)]
    if getattr(args, "te_txdb_rds", None):
        cmd += ["--txdb-rds", str(args.te_txdb_rds)]
    if getattr(args, "te_orgdb_package", None):
        cmd += ["--orgdb-package", str(args.te_orgdb_package)]

    step("Step 2/6 TE analysis: running atena/qtex + DESeq2 + ChIPseeker")
    run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))

    if not is_nonempty_file(final_csv):
        raise FileNotFoundError(f"Advanced TE analysis did not create expected final CSV: {final_csv}")
    return final_csv


def run_te_analysis_featurecounts(args, bam_samplesheet: str | Path, te_outdir: Path, final_csv: Path) -> Path:
    """Lightweight fallback TE mode: featureCounts + DESeq2 + simple TE-GTF metadata."""
    if not getattr(args, "te_gtf", None):
        raise ValueError("--te-mode simple requires --te-gtf. Advanced atena mode does not require this argument.")
    counts_dir = ensure_dir(te_outdir / "counts")
    de_dir = ensure_dir(te_outdir / "deseq2")
    annot_dir = ensure_dir(te_outdir / "annotation")

    df = pd.read_csv(bam_samplesheet, sep="\t")
    bams = df["bam_path"].astype(str).tolist()
    raw_counts = counts_dir / "TE_featureCounts.txt"
    log = Path(args.output_dir) / "pipeline_info" / "logs" / "TE_featureCounts.log"

    if not is_nonempty_file(raw_counts) or getattr(args, "force", False):
        cmd = [
            args.featurecounts_exe,
            "-T", str(args.threads),
            "-a", str(args.te_gtf),
            "-o", str(raw_counts),
            "-t", args.te_feature_type,
            "-g", args.te_attribute,
            "-s", _featurecounts_strand_code(args.strandedness),
        ]
        if args.paired:
            cmd += ["-p", "--countReadPairs"]
        cmd += bams
        step("Step 2/6 TE analysis: counting TE reads with featureCounts")
        run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))
    else:
        step(f"Step 2/6 TE analysis: reusing featureCounts table {raw_counts}")

    sample_meta = te_outdir / "samplesheet.bam.tsv"
    df.to_csv(sample_meta, sep="\t", index=False)
    rscript = Path(__file__).resolve().parents[1] / "r" / "te_deseq2_to_dsRNASeeker_csv.R"
    cmd = [
        args.rscript_exe, str(rscript),
        "--counts", str(raw_counts),
        "--samplesheet", str(sample_meta),
        "--case", args.case_label,
        "--control", args.control_label,
        "--te-gtf", str(args.te_gtf),
        "--out", str(final_csv),
        "--padj-max", str(args.te_padj_max),
        "--lfc-min", str(args.te_lfc_min),
    ]
    step("Step 2/6 TE analysis: running fallback featureCounts/DESeq2 annotation")
    run_cmd(cmd, log_path=Path(args.output_dir) / "pipeline_info" / "logs" / "TE_DESeq2_annotation.log", quiet=getattr(args, "quiet", True))
    return final_csv
