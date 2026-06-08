from __future__ import annotations

import json
import shutil
import types
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_dir, is_nonempty_file, step
from .alignment import build_star_index, run_star_alignment, markdup_bams, write_bam_samplesheet
from .te_analysis import run_te_analysis
from .rmats_runner import run_rmats_case_control
from .reditools_runner import run_reditools2
from .sprint_runner import run_sprint
from .run_pipeline import run_pipeline
from .summary import run_summary
from .delta import run_delta
from .zrna import run_zrna
from .inference import (
    infer_layout_and_read_length, infer_strandedness,
    apply_inferred_library_settings, infer_te_genome,
)


def _ns(**kwargs: Any) -> types.SimpleNamespace:
    return types.SimpleNamespace(**kwargs)


def _copy_public_args(args) -> dict[str, Any]:
    out = {}
    for k, v in vars(args).items():
        if isinstance(v, Path):
            out[k] = str(v)
        elif isinstance(v, (str, int, float, bool)) or v is None:
            out[k] = v
        elif isinstance(v, list):
            out[k] = [str(x) for x in v]
    return out


def _write_versions(args, outdir: Path) -> None:
    info = ensure_dir(outdir / "pipeline_info")
    checks = {
        "python": getattr(args, "python_exe", "python3"),
        "Rscript": getattr(args, "rscript_exe", "Rscript"),
        "STAR": getattr(args, "star_exe", "STAR"),
        "samtools": getattr(args, "samtools_exe", "samtools"),
        "bedtools": getattr(args, "bedtools_exe", "bedtools"),
        "featureCounts": getattr(args, "featurecounts_exe", "featureCounts"),
        "rmats": getattr(args, "rmats_exe", "rmats.py"),
        "RNAcofold": getattr(args, "rnacofold_exe", "RNAcofold"),
        "RNAfold": getattr(args, "rnafold_exe", "RNAfold"),
        "IntaRNA": getattr(args, "intarna_exe", "IntaRNA"),
    }
    lines = []
    for name, exe in checks.items():
        path = shutil.which(exe) or "NOT_FOUND"
        lines.append(f"{name}\t{exe}\t{path}")
    (info / "software_paths.tsv").write_text("tool\texecutable\tresolved_path\n" + "\n".join(lines) + "\n")
    (info / "dsRNASeeker_params.json").write_text(json.dumps(_copy_public_args(args), indent=2) + "\n")


def _normalize_fastq_samplesheet(path: str | Path, outdir: Path) -> pd.DataFrame:
    p = Path(path)
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep)
    required = {"sample_id", "condition"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"FASTQ samplesheet missing required columns: {sorted(missing)}")
    if "fastq_1" not in df.columns and "bam_path" not in df.columns:
        raise ValueError("Samplesheet must contain fastq_1 for FASTQ mode or bam_path for BAM mode.")
    if "fastq_2" not in df.columns:
        df["fastq_2"] = ""
    if "strandedness" not in df.columns:
        df["strandedness"] = "auto"
    df["sample_id"] = df["sample_id"].astype(str)
    df["condition"] = df["condition"].astype(str)
    return df


def _common_run_args(args, samplesheet_bam: Path, csv_in: Path, redit_dir: Path | None, sprint_dir: Path | None):
    return dict(
        output_dir=str(args.output_dir),
        case_label=args.case_label,
        control_label=args.control_label,
        samplesheet=str(samplesheet_bam),
        csv_in=str(csv_in),
        fasta=args.fasta,
        gtf=args.gtf,
        sprint_a2i_dir=str(sprint_dir) if sprint_dir else None,
        redit_dirs=[str(redit_dir)] if redit_dir else [],
        analyze_subset=args.analyze_subset,
        window_w=args.window_w,
        arm_aware=args.arm_aware,
        arm_pad=args.arm_pad,
        arm_min_cov=args.arm_min_cov,
        do_ddg=args.do_ddg,
        do_pf_interface=args.do_pf_interface,
        do_null_z=args.do_null_z,
        null_n=args.null_n,
        null_seed=args.null_seed,
        do_intarna=args.do_intarna,
        cofold_strong=args.cofold_strong,
        cofold_moderate=args.cofold_moderate,
        transcript_mapping_rscript=args.transcript_mapping_rscript,
        python_exe=args.python_exe,
        rscript_exe=args.rscript_exe,
        bedtools_exe=args.bedtools_exe,
        samtools_exe=args.samtools_exe,
        bamcoverage_exe=args.bamcoverage_exe,
        multibigwigsummary_exe=args.multibigwigsummary_exe,
        rnacofold_exe=args.rnacofold_exe,
        rnafold_exe=args.rnafold_exe,
        intarna_exe=args.intarna_exe,
        min_selected_candidates=args.min_selected_candidates,
    )


def run_workflow(args) -> None:
    outdir = ensure_dir(args.output_dir)
    args.output_dir = str(outdir)
    _write_versions(args, outdir)

    samples = _normalize_fastq_samplesheet(args.samplesheet, outdir)
    infer_layout_and_read_length(args, samples, outdir / 'pipeline_info' / 'inference')
    args.te_genome = infer_te_genome(args) if not args.skip_te_analysis and not args.precomputed_csv_in else getattr(args, 'te_genome', 'auto')
    case = args.case_label
    control = args.control_label
    observed = set(samples["condition"].astype(str))
    for label in [case, control]:
        if label not in observed:
            raise ValueError(f"condition label {label!r} not found in samplesheet conditions={sorted(observed)}")

    # Alignment / BAM contract
    if args.input_mode == "fastq":
        index_dir = Path(args.star_index) if args.star_index else outdir / "00_reference" / "star_index"
        # build_star_index() is resume-aware and validates STAR index sentinel files.
        build_star_index(args, index_dir)
        aligned = run_star_alignment(args, samples, index_dir, outdir / "01_alignment" / "star")
        bams = markdup_bams(args, aligned, outdir / "01_alignment" / "markdup")
        bam_sheet = write_bam_samplesheet(samples, bams, outdir / "pipeline_info" / "samplesheet.bam.tsv")
        inferred_strand = infer_strandedness(args, bam_sheet, outdir / 'pipeline_info' / 'inference')
        apply_inferred_library_settings(args, inferred_strand)
    else:
        if "bam_path" not in samples.columns:
            raise ValueError("--input-mode bam requires bam_path in samplesheet")
        bam_sheet = write_bam_samplesheet(samples, {r.sample_id: Path(r.bam_path) for r in samples.itertuples()}, outdir / "pipeline_info" / "samplesheet.bam.tsv")
        inferred_strand = infer_strandedness(args, bam_sheet, outdir / 'pipeline_info' / 'inference')
        apply_inferred_library_settings(args, inferred_strand)

    _write_versions(args, outdir)

    # TE CSV contract. --skip-te-analysis means "do not regenerate"; use an
    # explicit precomputed CSV first, otherwise reuse the expected internal CSV.
    expected_csv = outdir / "02_te" / "annotation" / f"TE_expression_annotation_{control}_vs_{case}_all_sig.dsRNASeeker.csv"
    if args.precomputed_csv_in:
        csv_in = Path(args.precomputed_csv_in)
        step(f"Step 2/6 TE analysis: using precomputed TE CSV {csv_in}")
    elif args.skip_te_analysis:
        if not is_nonempty_file(expected_csv):
            raise FileNotFoundError(f"--skip-te-analysis requested, but expected TE CSV is missing: {expected_csv}")
        csv_in = expected_csv
        step(f"Step 2/6 TE analysis: skipped; reusing {csv_in}")
    else:
        csv_in = run_te_analysis(args, bam_sheet, outdir / "02_te")

    # rMATS contract: b1=CASE, b2=CONTROL, so IncLevelDifference=CASE-CONTROL.
    expected_rmats_dir = outdir / "03_splicing" / "rmats"
    expected_rmats_file = expected_rmats_dir / f"RI.MATS.{args.rmats_track}.txt"
    if args.precomputed_rmats_dir:
        rmats_dir = Path(args.precomputed_rmats_dir)
        step(f"Step 3/6 splicing: using precomputed rMATS directory {rmats_dir}")
    elif args.skip_rmats:
        rmats_dir = expected_rmats_dir if is_nonempty_file(expected_rmats_file) else None
        step("Step 3/6 splicing: skipped" + (f"; reusing {rmats_dir}" if rmats_dir else "; no rMATS evidence will be added"))
    else:
        rmats_dir = run_rmats_case_control(args, bam_sheet, expected_rmats_dir)

    # RNA editing contracts. REDItools2 is default; SPRINT remains optional.
    expected_redit_dir = outdir / "04_editing" / "REDItools2"
    if args.precomputed_redit_dir:
        redit_dir = Path(args.precomputed_redit_dir)
        step(f"Step 4/6 RNA editing: using precomputed REDItools2 directory {redit_dir}")
    elif args.skip_reditools:
        redit_dir = expected_redit_dir if expected_redit_dir.exists() else None
        step("Step 4/6 RNA editing: REDItools2 skipped" + (f"; reusing {redit_dir}" if redit_dir else "; no REDItools2 evidence will be added"))
    else:
        redit_dir = run_reditools2(args, bam_sheet, expected_redit_dir)

    if args.run_sprint:
        if args.precomputed_sprint_dir:
            sprint_dir = Path(args.precomputed_sprint_dir)
            step(f"Step 4b/6 SPRINT: using precomputed SPRINT A-to-I directory {sprint_dir}")
        else:
            sprint_dir = run_sprint(args, samples, outdir / "04_editing" / "SPRINT")
    else:
        sprint_dir = Path(args.precomputed_sprint_dir) if args.precomputed_sprint_dir else None
        step("Step 4b/6 SPRINT: skipped by default")

    # Existing dsRNASeeker modules, now called internally.
    common = _common_run_args(args, bam_sheet, csv_in, redit_dir, sprint_dir)
    step("Step 5/6 dsRNASeeker core: running case condition")
    run_pipeline(_ns(**common, condition=case))
    step("Step 5/6 dsRNASeeker core: running control condition")
    run_pipeline(_ns(**common, condition=control))

    step("Step 6/6 summary/delta/Z-RNA: building fused summary")
    run_summary(_ns(
        output_dir=str(outdir), case_label=case, control_label=control, csv_in=str(csv_in),
        analyze_subset=args.analyze_subset, rmats_dir=str(rmats_dir) if rmats_dir else None,
        rmats_track=args.rmats_track, rmats_fdr_max=args.rmats_fdr_max,
        rmats_group1_label=case, rmats_group2_label=control, rmats_flip_dpsi=False,
        bedtools_exe=args.bedtools_exe, priority_top_n=args.priority_top_n,
        priority_mode=args.priority_mode, require_case_editing=args.require_case_editing,
        require_case_ri=args.require_case_ri, priority_score_mode=args.priority_score_mode,
        training_truth_table=args.training_truth_table, training_labels=args.training_labels,
        truth_symbol_col=args.truth_symbol_col, truth_label_mode=args.truth_label_mode,
        truth_label_col=args.truth_label_col, truth_padj_col=args.truth_padj_col,
        truth_logfc_col=args.truth_logfc_col, truth_padj_max=args.truth_padj_max,
        supervised_test_size=args.supervised_test_size, cv_folds=args.cv_folds,
        supervised_random_state=args.supervised_random_state,
    ))

    step("Step 6/6 summary/delta/Z-RNA: building delta table")
    run_delta(_ns(output_dir=str(outdir), case_label=case, control_label=control, analyze_subset=args.analyze_subset))
    step("Step 6/6 summary/delta/Z-RNA: annotating Z-RNA propensity")
    run_zrna(_ns(
        output_dir=str(outdir), case_label=case, control_label=control,
        analyze_subset=args.analyze_subset, summary_in=None, case_fasta=None, control_fasta=None,
        zrna_score_mode=args.zrna_score_mode, zrna_class_mode=args.zrna_class_mode,
        zrna_moderate_threshold=args.zrna_moderate_threshold, zrna_high_threshold=args.zrna_high_threshold,
    ))

    step(f"[OK] workflow complete: {outdir}")
