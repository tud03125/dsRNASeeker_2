from __future__ import annotations

from pathlib import Path
import shlex
import subprocess
import gzip
import shutil
import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


BWA_SUFFIXES = [".amb", ".ann", ".bwt", ".pac", ".sa"]


def _all_a2i_exist(samples: pd.DataFrame, a2i_dir: Path) -> bool:
    return all(is_nonempty_file(a2i_dir / f"{sid}_A_to_I.res") for sid in samples["sample_id"].astype(str))


def _bwa_index_files(fasta: str | Path) -> list[Path]:
    fasta = Path(fasta)
    return [Path(str(fasta) + suffix) for suffix in BWA_SUFFIXES]


def _bwa_index_complete(fasta: str | Path) -> bool:
    return all(is_nonempty_file(p) for p in _bwa_index_files(fasta))


def _materialize_sprint_fastq(path: str | Path, cache_dir: Path, auto_decompress: bool) -> str:
    src = Path(path)
    if not str(src).endswith((".gz", ".gzip")):
        return str(src)
    if not auto_decompress:
        raise ValueError(
            f"SPRINT does not accept gzipped FASTQ directly: {src}. "
            "Enable --sprint-auto-decompress or provide uncompressed FASTQ."
        )
    cache_dir.mkdir(parents=True, exist_ok=True)
    name = src.name
    for suffix in (".gzip", ".gz"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    dst = cache_dir / name
    if dst.exists() and dst.stat().st_size > 0 and dst.stat().st_mtime >= src.stat().st_mtime:
        return str(dst)
    step(f"Step 4b/6 SPRINT: decompressing {src.name} for legacy SPRINT input")
    tmp = Path(str(dst) + ".tmp")
    with gzip.open(src, "rb") as fin, tmp.open("wb") as fout:
        shutil.copyfileobj(fin, fout, length=16 * 1024 * 1024)
    tmp.replace(dst)
    return str(dst)


def ensure_sprint_bwa_index(args, fasta: str | Path, logs_dir: str | Path) -> None:
    """Build the BWA index needed by SPRINT if it is missing.

    SPRINT's FASTQ workflow uses BWA internally, so the reference FASTA needs the
    standard BWA index sidecar files: .amb, .ann, .bwt, .pac, .sa.
    """
    fasta = Path(fasta)
    logs_dir = ensure_dir(logs_dir)

    if _bwa_index_complete(fasta) and not getattr(args, "force", False):
        step(f"Step 4b/6 SPRINT/BWA index: reusing existing BWA index for {fasta}")
        return

    if getattr(args, "skip_sprint_bwa_index", False):
        missing = [str(p) for p in _bwa_index_files(fasta) if not is_nonempty_file(p)]
        if missing:
            raise FileNotFoundError(
                "--skip-sprint-bwa-index was used, but BWA index files are missing: " + ", ".join(missing)
            )
        return

    step(f"Step 4b/6 SPRINT/BWA index: building BWA index for {fasta}")
    cmd = [args.bwa_exe, "index", str(fasta)]
    run_cmd(cmd, log_path=logs_dir / "SPRINT_bwa_index.log", quiet=getattr(args, "quiet", True))


def ensure_sprint_prepare(args, outdir: str | Path) -> Path:
    """Run `sprint prepare` once per workflow unless already completed.

    SPRINT documents a `prepare` stage before `main`. The exact files created by
    SPRINT can vary by version, so dsRNASeeker records a local sentinel after the
    command succeeds. This keeps resume behavior stable without depending on
    SPRINT's internal file names.
    """
    outdir = ensure_dir(outdir)
    logs_dir = ensure_dir(Path(args.output_dir) / "pipeline_info" / "logs")
    sentinel = outdir / ".sprint_prepare.done"

    if sentinel.exists() and not getattr(args, "force", False):
        step(f"Step 4b/6 SPRINT prepare: reusing completed prepare step ({sentinel})")
        return sentinel

    if getattr(args, "skip_sprint_prepare", False):
        step("Step 4b/6 SPRINT prepare: skipped by user")
        return sentinel

    if not args.sprint_repeat_bed:
        raise ValueError("--run-sprint requires --sprint-repeat-bed; SPRINT uses repeat annotations via -rp.")

    ensure_sprint_bwa_index(args, args.fasta, logs_dir)

    cmd = [args.sprint_exe, "prepare", "-rp", str(args.sprint_repeat_bed), str(args.fasta), args.bwa_exe]
    if getattr(args, "sprint_prepare_extra", ""):
        cmd += shlex.split(args.sprint_prepare_extra)

    step("Step 4b/6 SPRINT prepare: running SPRINT reference/repeat preparation")
    run_cmd(cmd, log_path=logs_dir / "SPRINT_prepare.log", quiet=getattr(args, "quiet", True))
    sentinel.touch()
    return sentinel


def _run_geta2i(args, sample_out: Path, a2i: Path, sid: str) -> None:
    """Run getA2I.py or a wrapper around it.

    If --sprint-geta2i ends in .py, use --sprint-python-exe if provided, otherwise
    --python-exe. If it is a wrapper script (for example getA2I_py27.sh), execute
    the wrapper directly.
    """
    geta2i = str(args.sprint_geta2i)
    logs_dir = ensure_dir(Path(args.output_dir) / "pipeline_info" / "logs")

    if geta2i.endswith(".py"):
        py = getattr(args, "sprint_python_exe", None) or args.python_exe
        cmd = [py, geta2i, str(args.sprint_strand_specific), str(sample_out), str(a2i)]
    else:
        cmd = [geta2i, str(args.sprint_strand_specific), str(sample_out), str(a2i)]

    step(f"Step 4b/6 SPRINT: extracting A-to-I events for {sid}")
    run_cmd(cmd, log_path=logs_dir / f"SPRINT_getA2I.{sid}.log", quiet=getattr(args, "quiet", True))


def run_sprint(args, samples: pd.DataFrame, outdir: str | Path) -> Path:
    """Run generic two-stage SPRINT with explicit preparation/index support.

    Stage 0: BWA index + `sprint prepare` for the reference/repeat annotation.
    Stage 1: `sprint main` from FASTQ to sample-specific SPRINT output.
    Stage 2: getA2I extraction to <sample>_A_to_I.res for dsRNASeeker.

    This module is samplesheet-driven and contains no dataset-specific paths.
    """
    outdir = ensure_dir(outdir)
    raw_root = ensure_dir(outdir / "raw")
    a2i_dir = ensure_dir(outdir / "A_to_I")
    logs_dir = ensure_dir(Path(args.output_dir) / "pipeline_info" / "logs")

    if _all_a2i_exist(samples, a2i_dir) and not getattr(args, "force", False):
        step(f"Step 4b/6 SPRINT: reusing A-to-I outputs in {a2i_dir}")
        return a2i_dir

    if not args.sprint_repeat_bed:
        raise ValueError("--run-sprint requires --sprint-repeat-bed, because SPRINT main uses repeat annotations (-rp).")
    if not args.sprint_geta2i:
        raise ValueError("--run-sprint requires --sprint-geta2i pointing to SPRINT/utilities/getA2I.py or a wrapper script.")

    ensure_sprint_prepare(args, outdir / "prepare")
    fastq_cache = ensure_dir(outdir / "fastq_uncompressed")

    for row in samples.itertuples(index=False):
        sid = str(row.sample_id)
        sample_out = ensure_dir(raw_root / sid)
        done = sample_out / ".done"
        a2i = a2i_dir / f"{sid}_A_to_I.res"
        fq1 = _materialize_sprint_fastq(
            row.fastq_1, fastq_cache, getattr(args, 'sprint_auto_decompress', True)
        )
        raw_fq2 = str(getattr(row, "fastq_2", "") or "")
        fq2 = _materialize_sprint_fastq(
            raw_fq2, fastq_cache, getattr(args, 'sprint_auto_decompress', True)
        ) if raw_fq2 else ""

        if not done.exists() or getattr(args, "force", False):
            cmd = [
                args.sprint_exe, "main",
                "-rp", str(args.sprint_repeat_bed),
                str(args.fasta),
                "-1", fq1,
            ]
            if fq2:
                cmd += ["-2", fq2]
            cmd += [str(sample_out), args.bwa_exe, args.samtools_exe]
            if args.sprint_extra:
                cmd += shlex.split(args.sprint_extra)
            step(f"Step 4b/6 SPRINT: running main for {sid}")
            run_cmd(cmd, log_path=logs_dir / f"SPRINT_main.{sid}.log", quiet=getattr(args, "quiet", True))
            done.touch()
        else:
            step(f"Step 4b/6 SPRINT: reusing raw output for {sid}")

        if is_nonempty_file(a2i) and not getattr(args, "force", False):
            step(f"Step 4b/6 SPRINT: reusing A-to-I output for {sid}")
        else:
            _run_geta2i(args, sample_out, a2i, sid)

    return a2i_dir
