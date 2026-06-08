from __future__ import annotations

from pathlib import Path
import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _libtype_from_strandedness(s: str) -> str:
    s = (s or "auto").lower()
    if s in {"reverse", "fr-firststrand", "firststrand"}:
        return "fr-firststrand"
    if s in {"forward", "fr-secondstrand", "secondstrand"}:
        return "fr-secondstrand"
    return "fr-unstranded"


def run_rmats_case_control(args, bam_samplesheet: str | Path, rmats_outdir: str | Path) -> Path:
    rmats_outdir = ensure_dir(rmats_outdir)
    tmpdir = ensure_dir(Path(rmats_outdir) / "tmp")
    df = pd.read_csv(bam_samplesheet, sep="\t")
    case_bams = df.loc[df["condition"].astype(str) == args.case_label, "bam_path"].astype(str).tolist()
    ctrl_bams = df.loc[df["condition"].astype(str) == args.control_label, "bam_path"].astype(str).tolist()
    if not case_bams or not ctrl_bams:
        raise ValueError("rMATS requires at least one case and one control BAM")
    b1 = Path(rmats_outdir) / "b1_case.txt"
    b2 = Path(rmats_outdir) / "b2_control.txt"
    b1.write_text(",".join(case_bams) + "\n")
    b2.write_text(",".join(ctrl_bams) + "\n")
    out_file = Path(rmats_outdir) / f"RI.MATS.{args.rmats_track}.txt"
    if is_nonempty_file(out_file) and not getattr(args, "force", False):
        step(f"Step 3/6 splicing: reusing rMATS output {out_file}")
        return Path(rmats_outdir)
    libtype = args.rmats_libtype or _libtype_from_strandedness(args.strandedness)
    cmd = [
        args.rmats_exe,
        "--b1", str(b1),
        "--b2", str(b2),
        "--gtf", str(args.gtf),
        "-t", "paired" if args.paired else "single",
        "--libType", libtype,
        "--readLength", str(args.read_length),
        "--variable-read-length",
        "--allow-clipping",
        "--cstat", str(args.rmats_cstat),
        "--task", "both",
        "--nthread", str(args.threads),
        "--tstat", str(args.threads),
        "--od", str(rmats_outdir),
        "--tmp", str(tmpdir),
    ]
    if args.rmats_novel_ss:
        cmd.append("--novelSS")
    step("Step 3/6 splicing: running rMATS with b1=case and b2=control")
    run_cmd(cmd, log_path=Path(args.output_dir) / "pipeline_info" / "logs" / "rMATS.log", quiet=getattr(args, "quiet", True))
    return Path(rmats_outdir)
