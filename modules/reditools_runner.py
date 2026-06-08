from __future__ import annotations

from pathlib import Path
import pandas as pd
import shlex

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _all_filtered_exist(df: pd.DataFrame, outdir: Path) -> bool:
    return all(is_nonempty_file(outdir / f"{sid}_filtered_editing_events.txt") for sid in df["sample_id"].astype(str))


def run_reditools2(args, bam_samplesheet: str | Path, outdir: str | Path) -> Path:
    """Run a generic two-stage REDItools2 module.

    Stage 1: REDItools2/editing caller from BAM+FASTA -> raw per-sample table.
    Stage 2: generic A-to-I QC filter -> <sample>_filtered_editing_events.txt.

    No dataset-specific sample IDs or FCCC paths are hard-coded; all samples come
    from the BAM samplesheet created by the workflow.
    """
    outdir = ensure_dir(outdir)
    raw_dir = ensure_dir(outdir / "raw")
    filt_dir = ensure_dir(outdir)
    df = pd.read_csv(bam_samplesheet, sep="\t")

    if _all_filtered_exist(df, filt_dir) and not getattr(args, "force", False):
        step(f"Step 4/6 RNA editing: reusing REDItools2 filtered outputs in {filt_dir}")
        return filt_dir

    summaries = []
    for row in df.itertuples(index=False):
        sid = str(row.sample_id)
        cond = str(row.condition)
        bam = str(row.bam_path)
        raw = raw_dir / f"{sid}_reditools2_raw.txt"
        filt = filt_dir / f"{sid}_filtered_editing_events.txt"
        sample_summary = outdir / "qc" / f"{sid}_reditools2_summary.tsv"
        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"REDItools2.{sid}.log"

        if not is_nonempty_file(raw) or getattr(args, "force", False):
            # REDItools2 installations differ. This default matches the REDItools2
            # src/cineca/reditools.py convention used in your legacy scripts:
            #   reditools.py -f input.bam -r reference.fa -o output.txt -s <strand>
            cmd = [
                args.reditools_exe,
                "-f", bam,
                "-r", str(args.fasta),
                "-o", str(raw),
                "-s", str(args.reditools_strand),
            ]
            if args.reditools_extra:
                cmd += shlex.split(args.reditools_extra)
            step(f"Step 4/6 RNA editing: running REDItools2 for {sid}")
            run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))
        else:
            step(f"Step 4/6 RNA editing: reusing REDItools2 raw output for {sid}")

        if is_nonempty_file(filt) and not getattr(args, "force", False):
            step(f"Step 4/6 RNA editing: reusing REDItools2 filtered output for {sid}")
        else:
            rscript = Path(args.reditools_post_rscript) if args.reditools_post_rscript else Path(__file__).resolve().parents[1] / "r" / "reditools_filter_a2i.R"
            cmd = [
                args.rscript_exe, str(rscript),
                "--raw", str(raw),
                "--out", str(filt),
                "--sample", sid,
                "--condition", cond,
                "--strandedness", str(args.strandedness),
                "--min-meanq", str(args.reditools_min_meanq),
                "--min-coverage", str(args.reditools_min_coverage),
                "--min-frequency", str(args.reditools_min_frequency),
                "--summary-out", str(sample_summary),
            ]
            step(f"Step 4/6 RNA editing: filtering A-to-I REDItools2 events for {sid}")
            run_cmd(cmd, log_path=Path(args.output_dir) / "pipeline_info" / "logs" / f"REDItools2_filter.{sid}.log", quiet=getattr(args, "quiet", True))
        if sample_summary.exists():
            summaries.append(pd.read_csv(sample_summary, sep="\t"))

    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(outdir / "REDItools2_global_editing_index_summary.tsv", sep="\t", index=False)
    return filt_dir
