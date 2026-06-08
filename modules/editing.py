from __future__ import annotations
from pathlib import Path
import os
import re
import glob
import pandas as pd
from .utils import run_cmd


def _parse_reditools_sample_id(path: str) -> str:
    """
    Parse sample IDs from REDItools filtered_editing_events filenames.

    Handles:
      I1_filtered_editing_events.txt                  -> I1
      IA9_filtered_editing_events.txt                 -> IA9
      SRR1523658_filtered_editing_events.txt          -> SRR1523658
      2691149L1_filtered_editing_events.txt           -> 2691149L1
      thomagrp_..._stranded-1_2691149_filtered...txt  -> 2691149L1
      thomagrp_..._stranded-2_2691149_filtered...txt  -> 2691149L2
    """
    base = os.path.basename(path)

    # GSE308489 original REDItools filenames:
    # thomagrp_319096_RNAseq_total_stranded-1_2691149_filtered_editing_events.txt
    m = re.search(r"stranded-(?P<lane>[12])_(?P<num>\d+)_filtered_editing_events", base)
    if m:
        return f"{m.group('num')}L{m.group('lane')}"

    # Generic short-ID names:
    # I1_filtered_editing_events.txt
    # SRR1523658_filtered_editing_events.txt
    # 2691149L1_filtered_editing_events.txt
    base = re.sub(r"_filtered_editing_events.*$", "", base)

    return base


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """
    Case-insensitive column lookup while preserving the real column name.
    """
    exact = next((c for c in candidates if c in df.columns), None)
    if exact is not None:
        return exact

    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        hit = lower_map.get(c.lower())
        if hit is not None:
            return hit

    return None


def _write_event_bed_from_reditools(df: pd.DataFrame, out_bed: Path) -> int:
    """
    Convert REDItools event table to 0-based BED.

    Standard REDItools/HPC-REDItools columns are typically:
      Region, Position, Reference, Strand, ...
    where Region = chromosome and Position = 1-based genomic coordinate.
    """
    chr_col = _find_col(df, ["chrom", "chr", "Chromosome", "CHROM", "Region", "REGION"])
    pos_col = _find_col(df, ["position", "pos", "Position", "POS", "POSITION"])

    if chr_col is None or pos_col is None:
        raise ValueError(
            f"Could not find chromosome/position columns. "
            f"Columns={list(df.columns)}"
        )

    pos = pd.to_numeric(df[pos_col], errors="coerce")
    beddf = pd.DataFrame({
        "chrom": df[chr_col].astype(str),
        "start0": pos - 1,
        "end1": pos,
    })

    beddf = beddf.dropna(subset=["start0", "end1"])
    beddf["start0"] = beddf["start0"].astype(int)
    beddf["end1"] = beddf["end1"].astype(int)
    beddf = beddf[beddf["end1"] > 0].copy()

    if beddf.empty:
        return 0

    beddf[["chrom", "start0", "end1"]].to_csv(
        out_bed, sep="\t", header=False, index=False
    )
    return len(beddf)


def run_editing_overlays(args, outdir, tag, condition, cond_sample_ids):
    outdir = Path(outdir)
    bed = {
        "inverted": "pair_windows.inverted.bed",
        "hairpin": "pair_windows.hairpin.bed",
        "allpairs": "pair_windows.bed",
    }[tag]

    bed_path = outdir / bed
    cond_samples = set(str(x) for x in cond_sample_ids)

    if not bed_path.exists():
        print(f"[editing] WARNING: missing window BED: {bed_path}")
        return

    # -------------------------
    # SPRINT A-to-I overlays
    # -------------------------
    hit_rows = []
    if getattr(args, "sprint_a2i_dir", None) and Path(args.sprint_a2i_dir).is_dir():
        for f in glob.glob(os.path.join(args.sprint_a2i_dir, "*_A_to_I.res")):
            sid = os.path.basename(f).split("_A_to_I.res")[0]
            if cond_samples and sid not in cond_samples:
                continue

            try:
                df = pd.read_csv(
                    f,
                    sep="\t",
                    header=None,
                    names=["chrom", "start0", "end1", "type", "class", "strand", "ADDP"],
                    comment="#",
                )

                tmpbed = outdir / f".spr_{sid}.bed"
                df[["chrom", "start0", "end1"]].to_csv(
                    tmpbed, sep="\t", header=False, index=False
                )

                cp = run_cmd(
                    [
                        args.bedtools_exe,
                        "intersect",
                        "-c",
                        "-a",
                        str(bed_path),
                        "-b",
                        str(tmpbed),
                    ],
                    capture=True,
                )

                c = pd.read_csv(
                    pd.io.common.StringIO(cp.stdout),
                    sep="\t",
                    header=None,
                    names=["chrom", "start", "end", "pair_id", "cnt"],
                )
                c = c.groupby("pair_id")["cnt"].sum().reset_index()
                hit_rows.append(c)

            except Exception as e:
                print(f"[SPRINT] WARNING: failed on {f}: {e}")

    if hit_rows:
        H = pd.concat(hit_rows, ignore_index=True)
        W = H.groupby("pair_id")["cnt"].sum().reset_index()
        W.columns = ["pair_id", "AtoI_hits_window"]
        W.to_csv(outdir / f"AtoI_counts_window.{tag}.{condition}.tsv", sep="\t", index=False)
        print(f"[SPRINT] wrote AtoI_counts_window.{tag}.{condition}.tsv; total hits={int(W['AtoI_hits_window'].sum())}")
    else:
        print(f"[SPRINT] no A-to-I rows collected for condition={condition}, tag={tag}")

    # -------------------------
    # REDItools overlays
    # -------------------------
    hit_rows = []

    for d in getattr(args, "redit_dirs", None) or []:
        if not d or not Path(d).is_dir():
            print(f"[REDItools] WARNING: missing REDItools dir: {d}")
            continue

        files = (
            glob.glob(os.path.join(d, "*filtered_editing_events*.tsv"))
            + glob.glob(os.path.join(d, "*filtered_editing_events*.txt"))
        )

        print(f"[REDItools] scanning {d}; files={len(files)}; condition={condition}; allowed_samples={sorted(cond_samples)}")

        for f in files:
            sid = _parse_reditools_sample_id(f)

            if cond_samples and sid not in cond_samples:
                continue

            try:
                df = pd.read_csv(f, sep="\t")
                tmpbed = outdir / f".redi_{sid}.bed"

                n_events = _write_event_bed_from_reditools(df, tmpbed)
                if n_events == 0:
                    print(f"[REDItools] WARNING: zero valid positions after parsing {f}")
                    continue

                cp = run_cmd(
                    [
                        args.bedtools_exe,
                        "intersect",
                        "-c",
                        "-a",
                        str(bed_path),
                        "-b",
                        str(tmpbed),
                    ],
                    capture=True,
                )

                c = pd.read_csv(
                    pd.io.common.StringIO(cp.stdout),
                    sep="\t",
                    header=None,
                    names=["chrom", "start", "end", "pair_id", "cnt"],
                )

                c = c.groupby("pair_id")["cnt"].sum().reset_index()
                hit_rows.append(c)
                print(f"[REDItools] used sample={sid}; events={n_events}; file={os.path.basename(f)}")

            except Exception as e:
                print(f"[REDItools] WARNING: failed on {f}: {e}")

    if hit_rows:
        H = pd.concat(hit_rows, ignore_index=True)
        W = H.groupby("pair_id")["cnt"].sum().reset_index()
        W.columns = ["pair_id", "REDI_hits_window"]
        W.to_csv(outdir / f"REDI_counts_window.{tag}.{condition}.tsv", sep="\t", index=False)
        print(f"[REDItools] wrote REDI_counts_window.{tag}.{condition}.tsv; total hits={int(W['REDI_hits_window'].sum())}")
    else:
        print(f"[REDItools] no REDItools rows collected for condition={condition}, tag={tag}")
