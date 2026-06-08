from __future__ import annotations

from pathlib import Path
import pandas as pd

from .utils import ensure_dir, run_cmd, is_nonempty_file, step


def _star_index_ready(index_dir: Path) -> bool:
    # STAR creates these files in a usable genomeDir. Checking several avoids
    # treating a half-built/non-STAR directory as reusable.
    required = ["Genome", "SA", "SAindex", "genomeParameters.txt"]
    return index_dir.exists() and all((index_dir / x).exists() for x in required)


def build_star_index(args, index_dir: str | Path) -> Path:
    index_dir = ensure_dir(index_dir)
    if _star_index_ready(index_dir) and not getattr(args, "force", False):
        step(f"Step 1a/6 STAR index: reusing existing index at {index_dir}")
        return index_dir
    sjdb = int(args.sjdb_overhang) if args.sjdb_overhang is not None else max(int(args.read_length or 101) - 1, 1)
    log = Path(args.output_dir) / "pipeline_info" / "logs" / "STAR_genomeGenerate.log"
    cmd = [
        args.star_exe, "--runThreadN", str(args.threads),
        "--runMode", "genomeGenerate",
        "--genomeDir", str(index_dir),
        "--genomeFastaFiles", str(args.fasta),
        "--sjdbGTFfile", str(args.gtf),
        "--sjdbOverhang", str(sjdb),
    ]
    step(f"Step 1a/6 STAR index: building index at {index_dir}")
    run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))
    return index_dir


def run_star_alignment(args, samples: pd.DataFrame, index_dir: str | Path, star_outdir: str | Path) -> dict[str, Path]:
    star_outdir = ensure_dir(star_outdir)
    out: dict[str, Path] = {}
    for row in samples.itertuples(index=False):
        sid = str(row.sample_id)
        fq1 = str(row.fastq_1)
        fq2 = str(getattr(row, "fastq_2", "") or "")
        sample_prefix = star_outdir / f"{sid}_"
        bam = star_outdir / f"{sid}_Aligned.sortedByCoord.out.bam"
        if is_nonempty_file(bam) and not getattr(args, "force", False):
            step(f"Step 1b/6 STAR alignment: reusing {sid}")
            out[sid] = bam
            continue
        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"STAR_align.{sid}.log"
        cmd = [
            args.star_exe,
            "--runThreadN", str(args.threads),
            "--genomeDir", str(index_dir),
            "--readFilesIn", fq1,
        ]
        if fq2:
            cmd.append(fq2)
        if fq1.endswith((".gz", ".gzip")) or fq2.endswith((".gz", ".gzip")):
            cmd += ["--readFilesCommand", "zcat"]
        cmd += [
            "--outFileNamePrefix", str(sample_prefix),
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMstrandField", "intronMotif",
            "--quantMode", "GeneCounts",
        ]
        step(f"Step 1b/6 STAR alignment: running {sid}")
        run_cmd(cmd, log_path=log, quiet=getattr(args, "quiet", True))
        if not bam.exists():
            raise FileNotFoundError(f"STAR did not produce expected BAM: {bam}. See log: {log}")
        out[sid] = bam
    return out


def markdup_bams(args, aligned_bams: dict[str, Path], markdup_outdir: str | Path) -> dict[str, Path]:
    markdup_outdir = ensure_dir(markdup_outdir)
    out: dict[str, Path] = {}
    for sid, bam in aligned_bams.items():
        sorted_bam = markdup_outdir / f"{sid}.name_sorted.bam"
        fixmate_bam = markdup_outdir / f"{sid}.fixmate.bam"
        coord_bam = markdup_outdir / f"{sid}.coord_sorted.bam"
        md_bam = markdup_outdir / f"{sid}.markdup.sorted.bam"
        md_bai = Path(str(md_bam) + ".bai")
        if is_nonempty_file(md_bam) and is_nonempty_file(md_bai) and not getattr(args, "force", False):
            step(f"Step 1c/6 mark duplicates: reusing {sid}")
            out[sid] = md_bam
            continue
        log = Path(args.output_dir) / "pipeline_info" / "logs" / f"samtools_markdup.{sid}.log"
        step(f"Step 1c/6 mark duplicates: running {sid}")
        run_cmd([args.samtools_exe, "sort", "-@", str(args.threads), "-n", "-o", str(sorted_bam), str(bam)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "fixmate", "-m", str(sorted_bam), str(fixmate_bam)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "sort", "-@", str(args.threads), "-o", str(coord_bam), str(fixmate_bam)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "markdup", "-@", str(args.threads), str(coord_bam), str(md_bam)], log_path=log, quiet=getattr(args, "quiet", True))
        run_cmd([args.samtools_exe, "index", "-@", str(args.threads), str(md_bam)], log_path=log, quiet=getattr(args, "quiet", True))
        out[sid] = md_bam
    return out


def write_bam_samplesheet(samples: pd.DataFrame, bams: dict[str, Path], out_path: str | Path) -> Path:
    rows = []
    for row in samples.itertuples(index=False):
        sid = str(row.sample_id)
        if sid not in bams:
            raise ValueError(f"No BAM produced/found for sample {sid}")
        rows.append({"sample_id": sid, "condition": str(row.condition), "bam_path": str(bams[sid])})
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    return out_path
