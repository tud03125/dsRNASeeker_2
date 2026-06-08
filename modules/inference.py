from __future__ import annotations

import gzip
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

from .utils import ensure_dir, step


def _open_text(path: str | Path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith((".gz", ".gzip")) else open(p, "rt")


def _fastq_lengths(path: str | Path, records: int = 10000) -> list[int]:
    lengths: list[int] = []
    with _open_text(path) as fh:
        for _ in range(records):
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n\r")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n\r")
            if not plus or not qual:
                raise ValueError(f"Truncated FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ structure in {path}")
            if len(seq) != len(qual):
                raise ValueError(f"Sequence/quality length mismatch in {path}")
            lengths.append(len(seq))
    if not lengths:
        raise ValueError(f"No FASTQ records found in {path}")
    return lengths


def _bam_stats(path: str | Path, records: int = 10000) -> tuple[bool, list[int]]:
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("BAM inference requires pysam") from exc
    paired_votes = []
    lengths: list[int] = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.query_length is None:
                continue
            paired_votes.append(bool(read.is_paired))
            lengths.append(int(read.query_length))
            if len(lengths) >= records:
                break
    if not lengths:
        raise ValueError(f"No mapped reads available for inference in {path}")
    paired = sum(paired_votes) >= max(1, len(paired_votes) / 2)
    return paired, lengths


def infer_layout_and_read_length(args, samples: pd.DataFrame, info_dir: str | Path) -> None:
    """Infer paired/single layout and read length before STAR index generation."""
    info_dir = ensure_dir(info_dir)
    rows = []
    layout_votes: list[bool] = []
    all_lengths: list[int] = []
    inspect_n = int(getattr(args, "infer_fastq_records", 10000))

    if args.input_mode == "fastq":
        for row in samples.itertuples(index=False):
            fq1 = str(row.fastq_1)
            fq2 = str(getattr(row, "fastq_2", "") or "").strip()
            paired = bool(fq2)
            if not Path(fq1).exists():
                raise FileNotFoundError(fq1)
            if paired and not Path(fq2).exists():
                raise FileNotFoundError(fq2)
            lengths = _fastq_lengths(fq1, inspect_n)
            if paired:
                lengths += _fastq_lengths(fq2, inspect_n)
            layout_votes.append(paired)
            all_lengths.extend(lengths)
            rows.append({
                "sample_id": str(row.sample_id),
                "layout": "paired" if paired else "single",
                "min_read_length": min(lengths),
                "max_read_length": max(lengths),
                "modal_read_length": Counter(lengths).most_common(1)[0][0],
                "reads_inspected": len(lengths),
            })
    else:
        for row in samples.itertuples(index=False):
            paired, lengths = _bam_stats(row.bam_path, inspect_n)
            layout_votes.append(paired)
            all_lengths.extend(lengths)
            rows.append({
                "sample_id": str(row.sample_id),
                "layout": "paired" if paired else "single",
                "min_read_length": min(lengths),
                "max_read_length": max(lengths),
                "modal_read_length": Counter(lengths).most_common(1)[0][0],
                "reads_inspected": len(lengths),
            })

    if len(set(layout_votes)) != 1:
        raise ValueError("Mixed paired-end and single-end samples are not supported in one workflow run")
    inferred_paired = layout_votes[0]
    inferred_read_length = max(all_lengths)

    if getattr(args, "paired", None) is None:
        args.paired = inferred_paired
    elif bool(args.paired) != inferred_paired:
        raise ValueError(
            f"Explicit layout ({'paired' if args.paired else 'single'}) conflicts with input files "
            f"({'paired' if inferred_paired else 'single'})"
        )

    if getattr(args, "read_length", None) is None:
        args.read_length = inferred_read_length
    elif int(args.read_length) != inferred_read_length:
        step(
            f"[WARN] supplied --read-length {args.read_length} differs from observed maximum "
            f"{inferred_read_length}; retaining the explicit value"
        )

    pd.DataFrame(rows).to_csv(Path(info_dir) / "input_inference.tsv", sep="\t", index=False)
    step(
        f"Step 0/6 input inference: layout={'paired' if args.paired else 'single'}, "
        f"read_length={args.read_length}"
    )


def _parse_gtf_attrs(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, quoted, bare in re.findall(r"([A-Za-z0-9_.-]+)\s+(?:\"([^\"]*)\"|([^;\s]+))", text):
        out[key] = quoted or bare
    return out


def gtf_to_bed12(gtf: str | Path, bed: str | Path) -> Path:
    """Create transcript BED12 for RSeQC without external UCSC converters."""
    bed = Path(bed)
    if bed.exists() and bed.stat().st_size > 0:
        return bed
    transcripts: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    gene_names: dict[tuple[str, str, str], str] = {}
    with open(gtf, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2].lower() != "exon":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attrs = _parse_gtf_attrs(fields[8])
            tid = attrs.get("transcript_id") or attrs.get("transcriptId")
            if not tid:
                continue
            key = (chrom, strand, tid)
            transcripts[key].append((start - 1, end))
            gene_names[key] = attrs.get("gene_name") or attrs.get("gene_id") or tid
    if not transcripts:
        raise ValueError(f"No exon/transcript_id records could be parsed from GTF: {gtf}")
    bed.parent.mkdir(parents=True, exist_ok=True)
    with bed.open("wt") as out:
        for (chrom, strand, tid), exons in sorted(transcripts.items()):
            exons = sorted(set(exons))
            tx_start = min(x[0] for x in exons)
            tx_end = max(x[1] for x in exons)
            sizes = ",".join(str(e - s) for s, e in exons) + ","
            starts = ",".join(str(s - tx_start) for s, _ in exons) + ","
            name = f"{gene_names[(chrom, strand, tid)]}|{tid}"
            out.write(
                f"{chrom}\t{tx_start}\t{tx_end}\t{name}\t0\t{strand}\t"
                f"{tx_start}\t{tx_end}\t0\t{len(exons)}\t{sizes}\t{starts}\n"
            )
    return bed


def _parse_rseqc(text: str, paired: bool) -> tuple[float, float, float]:
    failed_match = re.search(r"Fraction of reads failed to determine:\s*([0-9.eE+-]+)", text)
    values = [float(x) for x in re.findall(r"Fraction of reads explained by .*?:\s*([0-9.eE+-]+)", text)]
    if not failed_match or len(values) < 2:
        raise ValueError(f"Could not parse infer_experiment.py output:\n{text}")
    return float(failed_match.group(1)), values[0], values[1]


def _classify_strand(first: float, second: float, stranded_threshold: float, unstranded_threshold: float) -> str:
    total = first + second
    if total <= 0:
        return "undetermined"
    f1, f2 = first / total, second / total
    if abs(f1 - f2) < unstranded_threshold:
        return "unstranded"
    if f1 >= stranded_threshold:
        return "forward"
    if f2 >= stranded_threshold:
        return "reverse"
    return "undetermined"


def infer_strandedness(args, bam_samplesheet: str | Path, info_dir: str | Path) -> str:
    """Infer library orientation from aligned BAMs using RSeQC."""
    explicit = str(getattr(args, "strandedness", "auto") or "auto").lower()
    if explicit != "auto":
        return explicit

    info_dir = ensure_dir(info_dir)
    bed = Path(getattr(args, "strandedness_bed", "") or (Path(info_dir) / "annotation.transcripts.bed12"))
    if not getattr(args, "strandedness_bed", None):
        gtf_to_bed12(args.gtf, bed)
    elif not bed.exists():
        raise FileNotFoundError(bed)

    df = pd.read_csv(bam_samplesheet, sep="\t")
    results = []
    executable = getattr(args, "infer_experiment_exe", "infer_experiment.py")
    sample_size = int(getattr(args, "strandedness_sample_size", 200000))
    stranded_threshold = float(getattr(args, "stranded_threshold", 0.8))
    unstranded_threshold = float(getattr(args, "unstranded_threshold", 0.1))

    for row in df.itertuples(index=False):
        cmd = [executable, "-r", str(bed), "-i", str(row.bam_path), "-s", str(sample_size)]
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = Path(info_dir) / f"infer_experiment.{row.sample_id}.log"
        log.write_text("# COMMAND\n" + " ".join(cmd) + "\n\n" + proc.stdout)
        if proc.returncode != 0:
            raise RuntimeError(f"RSeQC failed for {row.sample_id}; see {log}")
        failed, first, second = _parse_rseqc(proc.stdout, bool(args.paired))
        call = _classify_strand(first, second, stranded_threshold, unstranded_threshold)
        results.append({
            "sample_id": row.sample_id,
            "failed_fraction": failed,
            "orientation_1_fraction": first,
            "orientation_2_fraction": second,
            "inferred_strandedness": call,
        })

    calls = [r["inferred_strandedness"] for r in results]
    informative = [x for x in calls if x != "undetermined"]
    if informative and len(set(informative)) == 1:
        consensus = informative[0]
    elif informative:
        counts = Counter(informative)
        consensus, n = counts.most_common(1)[0]
        if n <= len(informative) / 2:
            consensus = "undetermined"
    else:
        consensus = "undetermined"

    if consensus == "undetermined":
        fallback = getattr(args, "strandedness_fallback", "unstranded")
        if fallback == "error":
            raise ValueError("Strandedness inference was undetermined; inspect pipeline_info/inference logs")
        step(f"[WARN] strandedness inference was undetermined; using fallback={fallback}")
        consensus = fallback

    pd.DataFrame(results).to_csv(Path(info_dir) / "strandedness_inference.tsv", sep="\t", index=False)
    (Path(info_dir) / "inferred_library.json").write_text(json.dumps({
        "paired": bool(args.paired),
        "read_length": int(args.read_length),
        "strandedness": consensus,
    }, indent=2) + "\n")
    step(f"Step 1d/6 strandedness inference: {consensus}")
    return consensus


def infer_te_genome(args) -> str:
    current = str(getattr(args, "te_genome", "auto") or "auto").lower()
    if current != "auto":
        return current
    hay = " ".join(str(x).lower() for x in [args.fasta, args.gtf, getattr(args, "te_rmsk_rds", "")])
    aliases = [
        ("hg38", ["hg38", "grch38"]),
        ("mm39", ["mm39", "grcm39"]),
        ("mm10", ["mm10", "grcm38"]),
    ]
    for value, tokens in aliases:
        if any(token in hay for token in tokens):
            step(f"Step 0/6 reference inference: te_genome={value}")
            return value
    if getattr(args, "te_txdb_gtf", None) and getattr(args, "te_rmsk_rds", None):
        step("Step 0/6 reference inference: te_genome=custom")
        return "custom"
    raise ValueError(
        "Could not infer --te-genome from reference filenames. Supply --te-genome "
        "or provide both --te-txdb-gtf and --te-rmsk-rds for a custom assembly."
    )


def apply_inferred_library_settings(args, strandedness: str) -> None:
    args.strandedness = strandedness
    if getattr(args, "rmats_libtype", None) is None:
        args.rmats_libtype = {
            "reverse": "fr-firststrand",
            "fr-firststrand": "fr-firststrand",
            "forward": "fr-secondstrand",
            "fr-secondstrand": "fr-secondstrand",
        }.get(strandedness, "fr-unstranded")
    if str(getattr(args, "reditools_strand", "auto")) == "auto":
        args.reditools_strand = {
            "reverse": "2",
            "fr-firststrand": "2",
            "forward": "1",
            "fr-secondstrand": "1",
        }.get(strandedness, "0")
    if str(getattr(args, "sprint_strand_specific", "auto")) == "auto":
        args.sprint_strand_specific = 0 if strandedness == "unstranded" else 1
    if getattr(args, "te_use_strand", None) is None:
        args.te_use_strand = strandedness != "unstranded"
