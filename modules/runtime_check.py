from __future__ import annotations
from pathlib import Path
from .utils import which_or_none, read_samplesheet


def check_runtime(args):
    msgs = []
    for exe_name in [args.python_exe, args.rscript_exe, args.bedtools_exe, args.samtools_exe,
                     args.bamcoverage_exe, args.multibigwigsummary_exe, args.rnacofold_exe,
                     args.rnafold_exe] + ([args.intarna_exe] if args.do_intarna else []):
        hit = which_or_none(exe_name)
        msgs.append(f"FOUND executable: {exe_name} -> {hit}" if hit else f"MISSING executable: {exe_name}")
    for p in [args.csv_in, args.fasta, args.gtf, args.samplesheet, args.transcript_mapping_rscript]:
        msgs.append(f"FOUND path: {p}" if p and Path(p).exists() else f"MISSING path: {p}")
    try:
        df = read_samplesheet(args.samplesheet)
        msgs.append(f"Samplesheet OK: {len(df)} rows; conditions={sorted(df['condition'].astype(str).unique().tolist())}")
    except Exception as e:
        msgs.append(f"Samplesheet ERROR: {e}")
    return msgs
